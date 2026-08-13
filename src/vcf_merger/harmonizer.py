"""Evidence-preserving VCF harmonization / ensemble aggregation."""

from __future__ import annotations

import heapq
import logging
import os
import tempfile
from collections import defaultdict
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pysam

from vcf_merger.callers import detect_caller, get_adapter, short_label
from vcf_merger.callers.base import CallerAdapter
from vcf_merger.canonical import is_supported_small_variant, make_canonical
from vcf_merger.exceptions import ValidationError
from vcf_merger.inspect import inspect_vcf
from vcf_merger.io_vcf import open_variant_file
from vcf_merger.models import (
    AnalysisMode,
    CallerEvidence,
    CanonicalVariant,
    EnsembleStrategy,
    HarmonizedVariant,
    InputVCFMetadata,
)
from vcf_merger.normalization import normalize_to_temp
from vcf_merger.provenance import build_provenance
from vcf_merger.validation import validate_input_set
from vcf_merger.writers import write_outputs

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class _SiteEvent:
    contig_idx: int
    pos: int
    ref: str
    alt: str
    source_idx: int
    evidence: CallerEvidence
    unsupported: bool = False
    unsupported_reason: Optional[str] = None

    def sort_key(self) -> tuple:
        return (self.contig_idx, self.pos, self.ref, self.alt, self.source_idx)


def _normalize_gt(gt: Optional[str]) -> Optional[str]:
    if gt is None:
        return None
    # Treat 0|1 and 0/1 as comparable by unphased allele multiset for conflict detection
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    try:
        alleles = tuple(sorted(int(p) if p != "." else -1 for p in parts))
    except ValueError:
        return gt
    return "/".join("." if a < 0 else str(a) for a in alleles)


def reconcile_genotypes(
    evidence: list[CallerEvidence],
    sample: str,
) -> tuple[Optional[str], bool]:
    """Return (representative_gt, conflict).

    Unanimous genotypes (ignoring phase) win; otherwise flag conflict and
    return a deterministic representative (first PASS, else first) genotype.
    """
    gts: list[tuple[str, bool, Optional[str]]] = []  # raw, passed, caller
    for e in evidence:
        for s in e.samples:
            if s.sample == sample and s.genotype:
                gts.append((s.genotype, e.passed, e.caller))
                break
    if not gts:
        return None, False
    norms = {_normalize_gt(g) for g, _, _ in gts}
    norms.discard(None)
    if len(norms) <= 1:
        # Prefer a PASS genotype's raw representation
        for raw, passed, _ in gts:
            if passed:
                return raw, False
        return gts[0][0], False
    for raw, passed, _ in gts:
        if passed:
            return raw, True
    return gts[0][0], True


def passes_strategy(
    evidence: list[CallerEvidence],
    strategy: EnsembleStrategy,
    *,
    consensus_n: int = 2,
    include_callers: Optional[set[str]] = None,
    exclude_callers: Optional[set[str]] = None,
) -> bool:
    filtered = evidence
    if include_callers is not None:
        filtered = [e for e in filtered if e.caller in include_callers]
    if exclude_callers:
        filtered = [e for e in filtered if e.caller not in exclude_callers]
    if not filtered:
        return False
    if strategy == EnsembleStrategy.UNION:
        return True
    if strategy == EnsembleStrategy.PASS_UNION:
        return any(e.passed for e in filtered)
    if strategy == EnsembleStrategy.CONSENSUS:
        return len({e.caller for e in filtered}) >= consensus_n
    if strategy == EnsembleStrategy.CALLER_SPECIFIC:
        return True
    return True


def _extract_from_record(
    record: pysam.VariantRecord,
    *,
    adapter: CallerAdapter,
    source_file: str,
    caller_version: Optional[str],
    detection_source,
    assembly: str,
    contig_index: dict[str, int],
    source_idx: int,
) -> Iterator[_SiteEvent]:
    alts = record.alts or ()
    if not alts:
        return
    for alt_index, alt in enumerate(alts):
        ok, reason = is_supported_small_variant(str(record.ref), str(alt))
        evidence = adapter.extract_evidence(
            record,
            alt_index=alt_index,
            source_file=source_file,
            caller_version=caller_version,
            detection_source=detection_source,
        )
        contig = str(record.contig)
        cidx = contig_index.get(contig, 10_000 + hash(contig) % 1000)
        yield _SiteEvent(
            contig_idx=cidx,
            pos=int(record.pos),
            ref=str(record.ref).upper(),
            alt=str(alt).upper(),
            source_idx=source_idx,
            evidence=evidence,
            unsupported=not ok,
            unsupported_reason=reason,
        )


def _stream_source_events(
    path: str,
    *,
    adapter: CallerAdapter,
    caller_version: Optional[str],
    detection_source,
    assembly: str,
    contig_index: dict[str, int],
    source_idx: int,
) -> Iterator[_SiteEvent]:
    """Yield events in merge-key order.

    Prefer contig-ordered ``fetch`` on indexed files; otherwise materialize and
    sort one caller stream (bounded by a single callset, not the cartesian product).
    """
    events: list[_SiteEvent] = []
    with open_variant_file(path) as vf:
        indexed = getattr(vf, "index", None) is not None
        if indexed and contig_index:
            # Iterate contigs in dictionary order for true streaming by contig.
            for contig, _idx in sorted(contig_index.items(), key=lambda x: x[1]):
                try:
                    iterator = vf.fetch(contig)
                except ValueError:
                    continue
                for record in iterator:
                    events.extend(
                        list(
                            _extract_from_record(
                                record,
                                adapter=adapter,
                                source_file=path,
                                caller_version=caller_version,
                                detection_source=detection_source,
                                assembly=assembly,
                                contig_index=contig_index,
                                source_idx=source_idx,
                            )
                        )
                    )
                # Flush per contig to keep memory bounded
                events.sort(key=lambda e: e.sort_key())
                yield from events
                events = []
            return

        for record in vf:
            events.extend(
                list(
                    _extract_from_record(
                        record,
                        adapter=adapter,
                        source_file=path,
                        caller_version=caller_version,
                        detection_source=detection_source,
                        assembly=assembly,
                        contig_index=contig_index,
                        source_idx=source_idx,
                    )
                )
            )
    events.sort(key=lambda e: e.sort_key())
    yield from events


def _merge_sorted_events(
    sources: list[Iterator[_SiteEvent]],
) -> Iterator[list[_SiteEvent]]:
    """Multi-way merge of sorted event streams; yield groups of identical alleles."""
    heap: list[tuple[tuple, int, _SiteEvent, Iterator[_SiteEvent]]] = []
    for i, it in enumerate(sources):
        try:
            ev = next(it)
        except StopIteration:
            continue
        heapq.heappush(heap, (ev.sort_key(), i, ev, it))

    current_key: Optional[tuple] = None
    group: list[_SiteEvent] = []

    while heap:
        key, i, ev, it = heapq.heappop(heap)
        allele_key = (ev.contig_idx, ev.pos, ev.ref, ev.alt)
        if current_key is None:
            current_key = allele_key
            group = [ev]
        elif allele_key == current_key:
            group.append(ev)
        else:
            yield group
            current_key = allele_key
            group = [ev]
        try:
            nxt = next(it)
            heapq.heappush(heap, (nxt.sort_key(), i, nxt, it))
        except StopIteration:
            continue
    if group:
        yield group


def harmonize_variant_group(
    events: list[_SiteEvent],
    *,
    assembly: str,
    strategy: EnsembleStrategy,
    consensus_n: int,
    include_callers: Optional[set[str]],
    exclude_callers: Optional[set[str]],
    primary_sample: str,
) -> Optional[HarmonizedVariant]:
    if not events:
        return None
    # Use first event locus
    e0 = events[0]
    chrom = e0.evidence.original_chrom or _contig_from_evidence(e0)
    canonical = make_canonical(assembly, chrom, e0.pos, e0.ref, e0.alt)
    evidence = [ev.evidence for ev in events]
    unsupported = any(ev.unsupported for ev in events)
    reason = next((ev.unsupported_reason for ev in events if ev.unsupported), None)

    if unsupported:
        return HarmonizedVariant(
            canonical=canonical,
            evidence=evidence,
            unsupported=True,
            unsupported_reason=reason,
        )

    if not passes_strategy(
        evidence,
        strategy,
        consensus_n=consensus_n,
        include_callers=include_callers,
        exclude_callers=exclude_callers,
    ):
        return None

    rep_gt, conflict = reconcile_genotypes(evidence, primary_sample)
    return HarmonizedVariant(
        canonical=canonical,
        evidence=evidence,
        gt_conflict=conflict,
        representative_gt=rep_gt,
    )


def _contig_from_evidence(ev: _SiteEvent) -> str:
    return ev.evidence.original_chrom or ""


def _build_contig_index(metas: list[InputVCFMetadata], reference: Optional[str]) -> dict[str, int]:
    """Contig order from reference dictionary when available, else first VCF header."""
    order: list[str] = []
    if reference and os.path.exists(reference + ".fai"):
        with open(reference + ".fai", encoding="utf-8") as fh:
            for line in fh:
                name = line.split("\t", 1)[0]
                if name not in order:
                    order.append(name)
    if not order:
        for c in metas[0].contigs:
            if c.name not in order:
                order.append(c.name)
    # include any extra contigs from other inputs
    for m in metas:
        for c in m.contigs:
            if c.name not in order:
                order.append(c.name)
    return {name: i for i, name in enumerate(order)}


def harmonize_vcfs(
    inputs: list[str | Path],
    output: str | Path,
    *,
    mode: AnalysisMode | str = AnalysisMode.GERMLINE,
    strategy: EnsembleStrategy | str = EnsembleStrategy.UNION,
    reference: Optional[str | Path] = None,
    normalize: bool = True,
    consensus_n: int = 2,
    include_callers: Optional[list[str]] = None,
    exclude_callers: Optional[list[str]] = None,
    callers: Optional[list[Optional[str]]] = None,
    tumor_sample: Optional[str] = None,
    normal_sample: Optional[str] = None,
    allow_gvcf: bool = False,
    command_line: Optional[list[str]] = None,
    keep_unsupported_in_sidecar: bool = True,
) -> dict[str, Any]:
    """Harmonize per-caller VCFs into VCF + evidence + provenance outputs."""
    mode = AnalysisMode(mode) if not isinstance(mode, AnalysisMode) else mode
    strategy = (
        EnsembleStrategy(strategy) if not isinstance(strategy, EnsembleStrategy) else strategy
    )
    input_paths = [str(p) for p in inputs]
    if len(input_paths) < 1:
        raise ValidationError("At least one input VCF is required")

    if callers is not None and len(callers) != len(input_paths):
        raise ValidationError("--caller count must match --input count when provided")

    # Inspect + detect callers
    metas: list[InputVCFMetadata] = []
    for i, path in enumerate(input_paths):
        explicit = callers[i] if callers else None
        metas.append(inspect_vcf(path, caller=explicit))

    ref_s = str(reference) if reference else None
    if normalize and not ref_s:
        raise ValidationError(
            "Normalization is enabled by default; provide --reference FASTA "
            "or pass normalize=False to skip (not recommended)."
        )

    assembly = validate_input_set(
        metas,
        mode=mode,
        reference=ref_s,
        allow_gvcf=allow_gvcf,
        require_identical_samples=True,
    )

    if mode == AnalysisMode.SOMATIC:
        _validate_somatic_samples(metas, tumor_sample=tumor_sample, normal_sample=normal_sample)

    # Normalize each input to temp when requested
    work_paths: list[str] = []
    temp_paths: list[str] = []
    try:
        if normalize and ref_s:
            for path in input_paths:
                tmp = normalize_to_temp(path, reference=ref_s)
                work_paths.append(tmp)
                temp_paths.append(tmp)
        else:
            work_paths = list(input_paths)

        contig_index = _build_contig_index(metas, ref_s)
        include_set = set(include_callers) if include_callers else None
        exclude_set = set(exclude_callers) if exclude_callers else None

        adapters: list[CallerAdapter] = []
        for m in metas:
            adapter = get_adapter(m.caller or "unknown")
            if adapter is None:
                # Unknown caller: use a generic FreeBayes-like passthrough via base subclass
                adapter = _GenericAdapter(m.caller or "unknown")
            adapters.append(adapter)

        sources = [
            _stream_source_events(
                work_paths[i],
                adapter=adapters[i],
                caller_version=metas[i].caller_version,
                detection_source=metas[i].caller_detection_source,
                assembly=assembly,
                contig_index=contig_index,
                source_idx=i,
            )
            for i in range(len(work_paths))
        ]

        samples = list(metas[0].samples)
        primary_sample = samples[0] if samples else "SAMPLE"
        if mode == AnalysisMode.SOMATIC and tumor_sample:
            primary_sample = tumor_sample

        variants: list[HarmonizedVariant] = []
        unsupported: list[HarmonizedVariant] = []
        for group in _merge_sorted_events(sources):
            # Fix contig name on canonical from evidence
            for ev in group:
                if not ev.evidence.original_chrom:
                    continue
            hv = harmonize_variant_group(
                group,
                assembly=assembly,
                strategy=strategy,
                consensus_n=consensus_n,
                include_callers=include_set,
                exclude_callers=exclude_set,
                primary_sample=primary_sample,
            )
            if hv is None:
                continue
            # Repair contig if empty
            if not hv.canonical.contig:
                chrom = group[0].evidence.original_chrom
                hv.canonical = make_canonical(
                    assembly, chrom, hv.canonical.pos, hv.canonical.ref, hv.canonical.alt
                )
            if hv.unsupported:
                unsupported.append(hv)
            else:
                variants.append(hv)

        # Collect headers from work paths
        input_headers: list[pysam.VariantHeader] = []
        for p in work_paths:
            with open_variant_file(p) as vf:
                input_headers.append(vf.header.copy())

        from vcf_merger import __version__

        tool_meta = [
            f"##vcf_mergerVersion={__version__}",
            f"##vcf_mergerCommand={' '.join(command_line or [])}",
            f"##vcf_mergerMode={mode.value}",
            f"##vcf_mergerStrategy={strategy.value}",
        ]
        if ref_s:
            tool_meta.append(f"##vcf_mergerReference={ref_s}")
        if mode == AnalysisMode.SOMATIC:
            if tumor_sample:
                tool_meta.append(f"##vcf_mergerTumorSample={tumor_sample}")
            if normal_sample:
                tool_meta.append(f"##vcf_mergerNormalSample={normal_sample}")

        prov = build_provenance(
            command_line=command_line or [],
            analysis_mode=mode.value,
            ensemble_strategy=strategy.value,
            reference_path=ref_s,
            normalization={
                "enabled": normalize,
                "tool": "bcftools norm" if normalize else None,
                "multiallelic": "-any" if normalize else None,
            },
            input_metas=metas,
            samples=samples,
        )

        all_for_sidecar = variants + (unsupported if keep_unsupported_in_sidecar else [])
        paths = write_outputs(
            output_vcf=output,
            variants=variants,
            input_headers=input_headers,
            samples=samples,
            provenance=prov,
            tool_meta_lines=tool_meta,
        )
        # Rewrite evidence sidecar to include unsupported if requested
        if keep_unsupported_in_sidecar and unsupported:
            from vcf_merger.writers import write_evidence_jsonl

            write_evidence_jsonl(paths["evidence"], all_for_sidecar)

        return {
            "paths": paths,
            "variant_count": len(variants),
            "unsupported_count": len(unsupported),
            "assembly": assembly,
            "samples": samples,
            "mode": mode.value,
            "strategy": strategy.value,
        }
    finally:
        for tmp in temp_paths:
            for suffix in ("", ".tbi", ".csi"):
                p = tmp + suffix
                if os.path.exists(p):
                    try:
                        os.unlink(p)
                    except OSError:
                        pass


def _validate_somatic_samples(
    metas: list[InputVCFMetadata],
    *,
    tumor_sample: Optional[str],
    normal_sample: Optional[str],
) -> None:
    all_samples = set()
    for m in metas:
        all_samples.update(m.samples)
    if tumor_sample and tumor_sample not in all_samples:
        raise ValidationError(
            f"--tumor-sample '{tumor_sample}' not found in inputs; available: {sorted(all_samples)}"
        )
    if normal_sample and normal_sample not in all_samples:
        raise ValidationError(
            f"--normal-sample '{normal_sample}' not found in inputs; available: {sorted(all_samples)}"
        )
    if not tumor_sample and len(all_samples) > 1:
        logger.warning(
            "Somatic mode without --tumor-sample; using first sample as representative."
        )


class _GenericAdapter(CallerAdapter):
    """Fallback adapter for unknown callers."""

    def __init__(self, name: str) -> None:
        self.name = name

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        return False, None
