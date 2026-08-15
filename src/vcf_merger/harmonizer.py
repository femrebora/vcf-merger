"""Evidence-preserving VCF harmonization / ensemble aggregation."""

from __future__ import annotations

import heapq
import logging
import os
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pysam

from vcf_merger.callers import get_adapter
from vcf_merger.callers.base import CallerAdapter
from vcf_merger.canonical import is_supported_small_variant, make_canonical
from vcf_merger.exceptions import ValidationError
from vcf_merger.inspect import inspect_vcf
from vcf_merger.io_vcf import open_variant_file
from vcf_merger.models import (
    AnalysisMode,
    CallerEvidence,
    EnsembleStrategy,
    HarmonizedVariant,
    InputVCFMetadata,
)
from vcf_merger.normalization import normalize_to_temp
from vcf_merger.provenance import build_provenance, write_provenance
from vcf_merger.somatic import SomaticContext, resolve_somatic_context
from vcf_merger.validation import validate_input_set
from vcf_merger.writers import StreamingHarmonizedWriter

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
    """Return (representative_gt, conflict)."""
    gts: list[tuple[str, bool, Optional[str]]] = []
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
    contig_index: dict[str, int],
    source_idx: int,
    tumor_sample: Optional[str] = None,
    normal_sample: Optional[str] = None,
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
            tumor_sample=tumor_sample,
            normal_sample=normal_sample,
        )
        contig = str(record.contig)
        cidx = contig_index.get(contig, 10_000 + abs(hash(contig)) % 1000)
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


def _is_indexed(path: str) -> bool:
    return (
        os.path.exists(path + ".tbi")
        or os.path.exists(path + ".csi")
        or path.endswith(".bcf")
    )


def _stream_source_events(
    path: str,
    *,
    adapter: CallerAdapter,
    caller_version: Optional[str],
    detection_source,
    contig_index: dict[str, int],
    source_idx: int,
    tumor_sample: Optional[str] = None,
    normal_sample: Optional[str] = None,
) -> Iterator[_SiteEvent]:
    """Yield events in merge-key order, flushing per contig when indexed."""

    def _emit(record: pysam.VariantRecord) -> list[_SiteEvent]:
        return list(
            _extract_from_record(
                record,
                adapter=adapter,
                source_file=path,
                caller_version=caller_version,
                detection_source=detection_source,
                contig_index=contig_index,
                source_idx=source_idx,
                tumor_sample=tumor_sample,
                normal_sample=normal_sample,
            )
        )

    with open_variant_file(path) as vf:
        if _is_indexed(path) and contig_index:
            for contig, _idx in sorted(contig_index.items(), key=lambda x: x[1]):
                events: list[_SiteEvent] = []
                try:
                    iterator = vf.fetch(contig)
                except ValueError:
                    continue
                for record in iterator:
                    events.extend(_emit(record))
                events.sort(key=lambda e: e.sort_key())
                yield from events
            return

        # Unindexed: one-caller materialization (still not a multi-caller cartesian join)
        events = []
        for record in vf:
            events.extend(_emit(record))
        events.sort(key=lambda e: e.sort_key())
        yield from events


def _merge_sorted_events(
    sources: list[Iterator[_SiteEvent]],
) -> Iterator[list[_SiteEvent]]:
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
        _key, i, ev, it = heapq.heappop(heap)
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
    e0 = events[0]
    chrom = e0.evidence.original_chrom
    canonical = make_canonical(assembly, chrom, e0.pos, e0.ref, e0.alt)
    evidence = [ev.evidence for ev in events]
    # Stable evidence order by caller name then source index
    evidence.sort(key=lambda e: (e.caller, e.source_file))
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


def _build_contig_index(metas: list[InputVCFMetadata], reference: Optional[str]) -> dict[str, int]:
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
    tumor_only: bool = False,
    allow_gvcf: bool = False,
    command_line: Optional[list[str]] = None,
    keep_unsupported_in_sidecar: bool = True,
    evidence_path: Optional[str | Path] = None,
    provenance_path: Optional[str | Path] = None,
) -> dict[str, Any]:
    """Harmonize per-caller VCFs into VCF + evidence + provenance outputs.

    Uses a multi-way merge over per-caller streams and writes incrementally so
    peak memory stays near one contig × N callers rather than full cartesian joins.
    """
    mode = AnalysisMode(mode) if not isinstance(mode, AnalysisMode) else mode
    strategy = (
        EnsembleStrategy(strategy) if not isinstance(strategy, EnsembleStrategy) else strategy
    )
    input_paths = [str(p) for p in inputs]
    if len(input_paths) < 1:
        raise ValidationError("At least one input VCF is required")

    if callers is not None and len(callers) != len(input_paths):
        raise ValidationError("--caller count must match --input count when provided")

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

    somatic_ctx: Optional[SomaticContext] = None
    if mode == AnalysisMode.SOMATIC:
        somatic_ctx = resolve_somatic_context(
            metas,
            tumor_sample=tumor_sample,
            normal_sample=normal_sample,
            tumor_only=tumor_only,
        )

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
                adapter = _GenericAdapter(m.caller or "unknown")
            adapters.append(adapter)

        tumor_name = somatic_ctx.tumor_sample if somatic_ctx else None
        normal_name = somatic_ctx.normal_sample if somatic_ctx else None

        sources = [
            _stream_source_events(
                work_paths[i],
                adapter=adapters[i],
                caller_version=metas[i].caller_version,
                detection_source=metas[i].caller_detection_source,
                contig_index=contig_index,
                source_idx=i,
                tumor_sample=tumor_name,
                normal_sample=normal_name,
            )
            for i in range(len(work_paths))
        ]

        samples = list(metas[0].samples)
        primary_sample = samples[0] if samples else "SAMPLE"
        if somatic_ctx is not None:
            primary_sample = somatic_ctx.tumor_sample

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
        if somatic_ctx is not None:
            tool_meta.append(f"##vcf_mergerSomatic={somatic_ctx.describe()}")
            tool_meta.append(f"##vcf_mergerTumorSample={somatic_ctx.tumor_sample}")
            if somatic_ctx.normal_sample:
                tool_meta.append(f"##vcf_mergerNormalSample={somatic_ctx.normal_sample}")

        out_base = str(output)
        if out_base.endswith(".vcf.gz"):
            stem = out_base[: -len(".vcf.gz")]
        elif out_base.endswith(".vcf"):
            stem = out_base[: -len(".vcf")]
        else:
            stem = out_base
        ev_path = str(evidence_path or f"{stem}.evidence.jsonl")
        prov_path = str(provenance_path or f"{stem}.provenance.json")

        with StreamingHarmonizedWriter(
            output,
            input_headers=input_headers,
            samples=samples,
            tool_meta_lines=tool_meta,
            evidence_path=ev_path,
            somatic=somatic_ctx,
        ) as writer:
            for group in _merge_sorted_events(sources):
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
                if hv.unsupported and not keep_unsupported_in_sidecar:
                    continue
                writer.write(hv, to_vcf=not hv.unsupported)

            paths = {
                "vcf": writer.output_vcf,
                "evidence": writer.evidence_path,
            }
            variant_count = writer.variant_count
            unsupported_count = writer.unsupported_count

        # close() already indexed; write provenance after stream completes
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
            notes=[
                "Technical VCF harmonization only; not ACMG/AMP classification.",
                "Caller concordance is not clinical evidence strength.",
            ]
            + (
                [
                    f"Somatic context: {somatic_ctx.describe()}",
                    "No AMP/ASCO/CAP clinical tier classification is performed.",
                ]
                if somatic_ctx
                else []
            ),
        )
        if somatic_ctx is not None:
            # Attach somatic roles into provenance inputs notes via normalization blob
            prov.normalization["somatic"] = {
                "tumor_sample": somatic_ctx.tumor_sample,
                "normal_sample": somatic_ctx.normal_sample,
                "tumor_only": somatic_ctx.tumor_only,
            }
        write_provenance(prov_path, prov)
        paths["provenance"] = prov_path

        return {
            "paths": paths,
            "variant_count": variant_count,
            "unsupported_count": unsupported_count,
            "assembly": assembly,
            "samples": samples,
            "mode": mode.value,
            "strategy": strategy.value,
            "somatic": (
                {
                    "tumor_sample": somatic_ctx.tumor_sample,
                    "normal_sample": somatic_ctx.normal_sample,
                    "tumor_only": somatic_ctx.tumor_only,
                }
                if somatic_ctx
                else None
            ),
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


class _GenericAdapter(CallerAdapter):
    """Fallback adapter for unknown callers."""

    def __init__(self, name: str) -> None:
        self.name = name

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        return False, None
