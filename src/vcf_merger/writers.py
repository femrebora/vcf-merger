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
from vcf_merger.somatic import SomaticContext

logger = logging.getLogger(__name__)


def evidence_to_dict(
    variant: HarmonizedVariant,
    *,
    somatic: Optional[SomaticContext] = None,
) -> dict[str, Any]:
    payload: dict[str, Any] = {
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
    if somatic is not None:
        payload["somatic"] = {
            "tumor_sample": somatic.tumor_sample,
            "normal_sample": somatic.normal_sample,
            "tumor_only": somatic.tumor_only,
        }
    return payload


def write_evidence_jsonl(
    path: str | Path,
    variants: Iterable[HarmonizedVariant],
    *,
    somatic: Optional[SomaticContext] = None,
) -> str:
    path_s = str(path)
    with open(path_s, "w", encoding="utf-8") as fh:
        for v in variants:
            fh.write(json.dumps(evidence_to_dict(v, somatic=somatic), sort_keys=True) + "\n")
    return path_s


def _representative_afs(
    variant: HarmonizedVariant,
) -> tuple[Optional[float], Optional[float]]:
    """Pick first non-null tumor/normal AF across evidence (PASS preferred)."""
    tumor = normal = None
    for prefer_pass in (True, False):
        for e in variant.evidence:
            if prefer_pass and not e.passed:
                continue
            if tumor is None and e.tumor_allele_fraction is not None:
                tumor = float(e.tumor_allele_fraction)
            if normal is None and e.normal_allele_fraction is not None:
                normal = float(e.normal_allele_fraction)
        if tumor is not None:
            break
    return tumor, normal


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


def _fill_record(
    out: pysam.VariantFile,
    v: HarmonizedVariant,
    samples: list[str],
    *,
    somatic: Optional[SomaticContext] = None,
) -> pysam.VariantRecord:
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
    if v.pass_caller_count > 0:
        rec.filter.add("PASS")
    else:
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

    if somatic is not None:
        pair = somatic.tumor_sample
        if somatic.normal_sample:
            pair = f"{somatic.tumor_sample}:{somatic.normal_sample}"
        rec.info["VM_SOMATIC_PAIR"] = pair
        taf, naf = _representative_afs(v)
        if taf is not None:
            rec.info["VM_TUMOR_AF"] = taf
        if naf is not None:
            rec.info["VM_NORMAL_AF"] = naf

    for sample in samples:
        gt = v.representative_gt if v.representative_gt else "./."
        dp = None
        gq = None
        for e in v.evidence:
            for s in e.samples:
                if s.sample == sample:
                    dp = s.depth if dp is None else dp
                    if s.genotype_quality is not None and gq is None:
                        gq = int(s.genotype_quality)
                    if v.representative_gt is None and s.genotype:
                        gt = s.genotype
                    break
        # Per-sample GT: prefer that sample's genotype when multi-sample
        sample_gt = None
        for e in v.evidence:
            for s in e.samples:
                if s.sample == sample and s.genotype:
                    sample_gt = s.genotype
                    break
            if sample_gt:
                break
        rec.samples[sample]["GT"] = _gt_tuple(sample_gt or gt)
        if dp is not None:
            rec.samples[sample]["DP"] = dp
        if gq is not None:
            rec.samples[sample]["GQ"] = gq
    return rec


class StreamingHarmonizedWriter:
    """Write VCF + evidence JSONL incrementally (WGS-friendly)."""

    def __init__(
        self,
        output_vcf: str | Path,
        *,
        input_headers: list[pysam.VariantHeader],
        samples: list[str],
        tool_meta_lines: Optional[list[str]] = None,
        evidence_path: Optional[str | Path] = None,
        somatic: Optional[SomaticContext] = None,
    ) -> None:
        self.output_vcf = str(output_vcf)
        self.samples = samples
        self.somatic = somatic
        base = self.output_vcf
        if base.endswith(".vcf.gz"):
            stem = base[: -len(".vcf.gz")]
        elif base.endswith(".vcf"):
            stem = base[: -len(".vcf")]
        else:
            stem = base
        self.evidence_path = str(evidence_path or f"{stem}.evidence.jsonl")
        self.header = reconcile_headers(
            input_headers,
            samples=samples,
            tool_meta_lines=tool_meta_lines,
        )
        mode = "wz" if self.output_vcf.endswith(".gz") else "w"
        parent = os.path.dirname(os.path.abspath(self.output_vcf))
        if parent:
            os.makedirs(parent, exist_ok=True)
        self._vcf = pysam.VariantFile(self.output_vcf, mode, header=self.header)
        self._evidence_fh = open(self.evidence_path, "w", encoding="utf-8")
        self.variant_count = 0
        self.unsupported_count = 0

    def write(self, variant: HarmonizedVariant, *, to_vcf: bool = True) -> None:
        if variant.unsupported:
            self.unsupported_count += 1
            self._evidence_fh.write(
                json.dumps(evidence_to_dict(variant, somatic=self.somatic), sort_keys=True) + "\n"
            )
            return
        if to_vcf:
            rec = _fill_record(self._vcf, variant, self.samples, somatic=self.somatic)
            self._vcf.write(rec)
            self.variant_count += 1
        self._evidence_fh.write(
            json.dumps(evidence_to_dict(variant, somatic=self.somatic), sort_keys=True) + "\n"
        )

    def close(self) -> dict[str, str]:
        self._vcf.close()
        self._evidence_fh.close()
        if self.output_vcf.endswith(".gz"):
            try:
                index_variant_file(self.output_vcf)
            except Exception as exc:  # noqa: BLE001
                logger.warning("Indexing failed for %s: %s", self.output_vcf, exc)
        return {"vcf": self.output_vcf, "evidence": self.evidence_path}

    def __enter__(self) -> "StreamingHarmonizedWriter":
        return self

    def __exit__(self, exc_type, exc, tb) -> None:
        self.close()


def write_harmonized_vcf(
    path: str | Path,
    variants: list[HarmonizedVariant],
    *,
    input_headers: list[pysam.VariantHeader],
    samples: list[str],
    tool_meta_lines: Optional[list[str]] = None,
    write_unsupported: bool = False,
    somatic: Optional[SomaticContext] = None,
) -> str:
    """Write a compact standards-compliant VCF/VCF.GZ (batch helper)."""
    with StreamingHarmonizedWriter(
        path,
        input_headers=input_headers,
        samples=samples,
        tool_meta_lines=tool_meta_lines,
        somatic=somatic,
    ) as writer:
        for v in variants:
            if v.unsupported and not write_unsupported:
                writer.write(v, to_vcf=False)
            else:
                writer.write(v, to_vcf=not v.unsupported)
        return writer.output_vcf


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
    somatic: Optional[SomaticContext] = None,
) -> dict[str, str]:
    base = str(output_vcf)
    if base.endswith(".vcf.gz"):
        stem = base[: -len(".vcf.gz")]
    elif base.endswith(".vcf"):
        stem = base[: -len(".vcf")]
    else:
        stem = base
    ev_path = str(evidence_path or f"{stem}.evidence.jsonl")
    prov_path = str(provenance_path or f"{stem}.provenance.json")

    with StreamingHarmonizedWriter(
        output_vcf,
        input_headers=input_headers,
        samples=samples,
        tool_meta_lines=tool_meta_lines,
        evidence_path=ev_path,
        somatic=somatic,
    ) as writer:
        for v in variants:
            writer.write(v, to_vcf=not v.unsupported)

    write_provenance(prov_path, provenance)
    return {"vcf": str(output_vcf), "evidence": ev_path, "provenance": prov_path}
