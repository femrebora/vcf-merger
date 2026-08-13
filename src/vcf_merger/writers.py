"""Output writers for harmonized VCF, evidence sidecar, and provenance."""

from __future__ import annotations

import json
import logging
import os
from pathlib import Path
from typing import Any, Iterable, Optional

import pysam

from vcf_merger.callers import short_label
from vcf_merger.headers import reconcile_headers
from vcf_merger.io_vcf import index_variant_file
from vcf_merger.models import HarmonizedVariant, ProvenanceRecord
from vcf_merger.provenance import write_provenance

logger = logging.getLogger(__name__)


def evidence_to_dict(variant: HarmonizedVariant) -> dict[str, Any]:
    return {
        "canonical": {
            "assembly": variant.canonical.assembly,
            "contig": variant.canonical.contig,
            "pos": variant.canonical.pos,
            "ref": variant.canonical.ref,
            "alt": variant.canonical.alt,
            "vrs_id": variant.canonical.vrs_id,
        },
        "callers_all": variant.callers_all,
        "callers_pass": variant.callers_pass,
        "caller_count": variant.caller_count,
        "pass_caller_count": variant.pass_caller_count,
        "gt_conflict": variant.gt_conflict,
        "representative_gt": variant.representative_gt,
        "unsupported": variant.unsupported,
        "unsupported_reason": variant.unsupported_reason,
        "evidence": [
            {
                "caller": e.caller,
                "caller_version": e.caller_version,
                "source_file": e.source_file,
                "filter_status": e.filter_status,
                "passed": e.passed,
                "qual": e.qual,
                "original": {
                    "chrom": e.original_chrom,
                    "pos": e.original_pos,
                    "ref": e.original_ref,
                    "alt": e.original_alt,
                },
                "caller_specific_metrics": e.caller_specific_metrics,
                "raw_info": e.raw_info,
                "tumor_allele_fraction": e.tumor_allele_fraction,
                "normal_allele_fraction": e.normal_allele_fraction,
                "samples": [
                    {
                        "sample": s.sample,
                        "genotype": s.genotype,
                        "phased": s.phased,
                        "ploidy": s.ploidy,
                        "depth": s.depth,
                        "ref_depth": s.ref_depth,
                        "alt_depth": s.alt_depth,
                        "allele_fraction": s.allele_fraction,
                        "genotype_quality": s.genotype_quality,
                        "raw_format": s.raw_format,
                    }
                    for s in e.samples
                ],
            }
            for e in variant.evidence
        ],
    }


def write_evidence_jsonl(path: str | Path, variants: Iterable[HarmonizedVariant]) -> str:
    path_s = str(path)
    with open(path_s, "w", encoding="utf-8") as fh:
        for v in variants:
            fh.write(json.dumps(evidence_to_dict(v), sort_keys=True) + "\n")
    return path_s


def write_harmonized_vcf(
    path: str | Path,
    variants: list[HarmonizedVariant],
    *,
    input_headers: list[pysam.VariantHeader],
    samples: list[str],
    tool_meta_lines: Optional[list[str]] = None,
    write_unsupported: bool = False,
) -> str:
    """Write a compact standards-compliant VCF/VCF.GZ."""
    path_s = str(path)
    header = reconcile_headers(
        input_headers,
        samples=samples,
        tool_meta_lines=tool_meta_lines,
    )

    mode = "wz" if path_s.endswith(".gz") else "w"
    # Sort variants by contig order in header then position
    contig_order = {name: i for i, name in enumerate(header.contigs)}

    def sort_key(v: HarmonizedVariant) -> tuple[int, int, str, str]:
        return (
            contig_order.get(v.canonical.contig, 10_000),
            v.canonical.pos,
            v.canonical.ref,
            v.canonical.alt,
        )

    ordered = sorted(
        (v for v in variants if write_unsupported or not v.unsupported),
        key=sort_key,
    )

    with pysam.VariantFile(path_s, mode, header=header) as out:
        for v in ordered:
            rec = out.new_record(
                contig=v.canonical.contig,
                start=v.canonical.pos - 1,
                stop=v.canonical.pos - 1 + len(v.canonical.ref),
                alleles=(v.canonical.ref, v.canonical.alt),
                id=".",
                qual=None,
                filter=None,
                info={},
            )
            # FILTER: PASS if any pass caller else set to first non-pass or '.'
            if v.pass_caller_count > 0:
                rec.filter.add("PASS")
            else:
                # leave unfiltered / use first filter token if available
                statuses = [e.filter_status for e in v.evidence if e.filter_status not in (".", "PASS")]
                if statuses:
                    for tok in statuses[0].split(";"):
                        if tok and tok not in rec.filter:
                            try:
                                rec.filter.add(tok)
                            except Exception:  # noqa: BLE001
                                pass

            labels_all = [short_label(c) for c in v.callers_all]
            labels_pass = [short_label(c) for c in v.callers_pass]
            rec.info["VM_CALLERS"] = ",".join(labels_all)
            rec.info["VM_PASS_CALLERS"] = ",".join(labels_pass) if labels_pass else "."
            rec.info["VM_CALLER_COUNT"] = v.caller_count
            rec.info["VM_PASS_CALLER_COUNT"] = v.pass_caller_count
            if v.gt_conflict:
                rec.info["VM_GT_CONFLICT"] = True
            if v.evidence:
                e0 = v.evidence[0]
                rec.info["VM_ORIG"] = (
                    f"{e0.original_chrom}:{e0.original_pos}:{e0.original_ref}>{e0.original_alt}"
                )

            # Representative sample FORMAT
            # Prefer first sample in samples list
            for sample in samples:
                gt = v.representative_gt if not v.gt_conflict else (
                    v.representative_gt if v.representative_gt else "./."
                )
                # Gather DP/GQ from first evidence that has this sample
                dp = None
                gq = None
                for e in v.evidence:
                    for s in e.samples:
                        if s.sample == sample:
                            if gt is None or (v.gt_conflict and s.genotype):
                                # keep representative chosen earlier
                                pass
                            dp = s.depth if dp is None else dp
                            gq = int(s.genotype_quality) if s.genotype_quality is not None and gq is None else gq
                            if v.representative_gt is None and s.genotype:
                                gt = s.genotype
                            break
                rec.samples[sample]["GT"] = _gt_tuple(gt)
                if dp is not None:
                    rec.samples[sample]["DP"] = dp
                if gq is not None:
                    rec.samples[sample]["GQ"] = gq

            out.write(rec)

    if path_s.endswith(".gz"):
        try:
            index_variant_file(path_s)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Indexing failed for %s: %s", path_s, exc)

    return path_s


def _gt_tuple(gt: Optional[str]) -> tuple:
    if not gt or gt == ".":
        return (None,)
    sep = "|" if "|" in gt else "/"
    parts = gt.split(sep)
    out = []
    for p in parts:
        if p == ".":
            out.append(None)
        else:
            try:
                out.append(int(p))
            except ValueError:
                out.append(None)
    return tuple(out)


def write_outputs(
    *,
    output_vcf: str | Path,
    variants: list[HarmonizedVariant],
    input_headers: list[pysam.VariantHeader],
    samples: list[str],
    provenance: ProvenanceRecord,
    tool_meta_lines: Optional[list[str]] = None,
    evidence_path: Optional[str | Path] = None,
    provenance_path: Optional[str | Path] = None,
) -> dict[str, str]:
    vcf_path = write_harmonized_vcf(
        output_vcf,
        variants,
        input_headers=input_headers,
        samples=samples,
        tool_meta_lines=tool_meta_lines,
    )
    base = str(output_vcf)
    if base.endswith(".vcf.gz"):
        stem = base[: -len(".vcf.gz")]
    elif base.endswith(".vcf"):
        stem = base[: -len(".vcf")]
    else:
        stem = base

    ev_path = str(evidence_path or f"{stem}.evidence.jsonl")
    prov_path = str(provenance_path or f"{stem}.provenance.json")
    write_evidence_jsonl(ev_path, variants)
    write_provenance(prov_path, provenance)
    return {"vcf": vcf_path, "evidence": ev_path, "provenance": prov_path}
