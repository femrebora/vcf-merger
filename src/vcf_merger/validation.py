"""Input validation: assembly, contigs, samples, gVCF, mode compatibility."""

from __future__ import annotations

import logging
import os
import re
from pathlib import Path
from typing import Optional

import pysam

from vcf_merger.exceptions import (
    GVcfRejectedError,
    ReferenceError,
    SampleMismatchError,
    ValidationError,
)
from vcf_merger.io_vcf import (
    contig_list_from_header,
    detect_file_kind,
    header_assembly_hints,
    open_variant_file,
)
from vcf_merger.models import AnalysisMode, ContigInfo, FileKind, InputVCFMetadata

logger = logging.getLogger(__name__)

_NON_REF = re.compile(r"<NON_REF>", re.IGNORECASE)


def detect_gvcf_signals(
    path: str | Path,
    header: pysam.VariantHeader,
    *,
    max_scan_records: int = 200,
) -> tuple[bool, list[str]]:
    """Detect gVCF / reference-confidence characteristics."""
    reasons: list[str] = []
    kind = detect_file_kind(path)
    if kind in (FileKind.GVCF, FileKind.GVCF_GZ):
        reasons.append("filename_suggests_gvcf")

    # ALT header / ##ALT=<ID=NON_REF
    header_text = str(header)
    if "NON_REF" in header_text:
        reasons.append("header_contains_NON_REF")
    if 'ID=END' in header_text or "INFO=<ID=END" in header_text:
        # END alone is not definitive but with NON_REF is strong
        if "NON_REF" in header_text:
            reasons.append("header_END_with_NON_REF")

    # Scan early records
    try:
        with open_variant_file(path) as vf:
            for i, rec in enumerate(vf):
                if i >= max_scan_records:
                    break
                alts = rec.alts or ()
                if any(a == "<NON_REF>" or a == "<*>" for a in alts):
                    reasons.append("record_has_NON_REF_or_symbolic_ref")
                    break
                if "END" in rec.info:
                    # reference block style
                    if not alts or all(a in (".", "<NON_REF>", "<*>") for a in alts):
                        reasons.append("reference_confidence_block_END")
                        break
    except Exception as exc:  # noqa: BLE001
        logger.debug("gVCF scan incomplete for %s: %s", path, exc)

    # Deduplicate while preserving order
    seen: set[str] = set()
    uniq: list[str] = []
    for r in reasons:
        if r not in seen:
            seen.add(r)
            uniq.append(r)

    is_gvcf = bool(uniq) and (
        "filename_suggests_gvcf" in uniq
        or "record_has_NON_REF_or_symbolic_ref" in uniq
        or "header_contains_NON_REF" in uniq
        or "reference_confidence_block_END" in uniq
    )
    return is_gvcf, uniq


def collect_metadata(path: str | Path) -> InputVCFMetadata:
    """Inspect a VCF/BCF and return structured metadata (no full load)."""
    path_s = str(path)
    kind = detect_file_kind(path_s)
    with open_variant_file(path_s) as vf:
        header = vf.header
        assembly, reference = header_assembly_hints(header)
        contigs = contig_list_from_header(header)
        samples = list(header.samples)
        info_ids = list(header.info.keys())
        format_ids = list(header.formats.keys())
        filter_ids = list(header.filters.keys())
        fileformat = None
        for rec in header.records:
            if getattr(rec, "key", "") == "fileformat":
                fileformat = getattr(rec, "value", None)
                break
        is_gvcf, reasons = detect_gvcf_signals(path_s, header)

    return InputVCFMetadata(
        path=path_s,
        file_kind=kind,
        fileformat=fileformat,
        assembly=assembly,
        reference=reference,
        contigs=contigs,
        samples=samples,
        is_gvcf=is_gvcf,
        gvcf_reasons=reasons,
        header_info_ids=info_ids,
        header_format_ids=format_ids,
        header_filter_ids=filter_ids,
    )


def reject_gvcf_if_needed(meta: InputVCFMetadata, *, allow_gvcf: bool = False) -> None:
    if meta.is_gvcf and not allow_gvcf:
        reasons = ", ".join(meta.gvcf_reasons) or "gVCF characteristics detected"
        raise GVcfRejectedError(
            f"Refusing to ensemble-merge gVCF/reference-confidence file '{meta.path}' "
            f"({reasons}). GATK HaplotypeCaller gVCFs are intermediate representations "
            f"and should be genotyped (e.g. GenotypeGVCFs) before variant-level "
            f"ensemble merging. Pass an explicit allow flag only if a future gVCF "
            f"module is intentionally enabled."
        )


def validate_reference_fasta(reference: str | Path) -> str:
    path = str(reference)
    if not os.path.exists(path):
        raise ReferenceError(f"Reference FASTA not found: {path}")
    # Accept plain or bgzipped fasta; require .fai if present or buildable
    fai = path + ".fai"
    if not os.path.exists(fai):
        try:
            pysam.faidx(path)
        except Exception as exc:  # noqa: BLE001
            raise ReferenceError(
                f"Reference FASTA index missing and could not be created for '{path}': {exc}"
            ) from exc
    return path


def infer_assembly_from_reference(reference: str | Path) -> Optional[str]:
    name = os.path.basename(str(reference)).lower()
    if "grch38" in name or "hg38" in name or "hs38" in name:
        return "GRCh38"
    if "grch37" in name or "hg19" in name or "hs37" in name:
        return "GRCh37"
    return None


def contig_dict_compatible(
    left: list[ContigInfo],
    right: list[ContigInfo],
) -> tuple[bool, list[str]]:
    """Check contig name/length compatibility for shared contigs."""
    problems: list[str] = []
    right_map = {c.name: c for c in right}
    for c in left:
        if c.name not in right_map:
            continue
        other = right_map[c.name]
        if c.length is not None and other.length is not None and c.length != other.length:
            problems.append(
                f"contig length mismatch for {c.name}: {c.length} vs {other.length}"
            )
    return (len(problems) == 0), problems


def validate_input_set(
    metas: list[InputVCFMetadata],
    *,
    mode: AnalysisMode,
    reference: Optional[str] = None,
    allow_gvcf: bool = False,
    require_identical_samples: bool = True,
) -> str:
    """Validate a set of inputs before merge. Returns resolved assembly label."""
    if not metas:
        raise ValidationError("No input VCFs provided")

    for m in metas:
        reject_gvcf_if_needed(m, allow_gvcf=allow_gvcf)

    # Assembly consistency
    assemblies = {m.assembly for m in metas if m.assembly}
    ref_assembly = infer_assembly_from_reference(reference) if reference else None
    if reference:
        validate_reference_fasta(reference)

    if len(assemblies) > 1:
        raise ValidationError(
            f"Incompatible reference assemblies across inputs: {sorted(assemblies)}. "
            f"Never silently merge GRCh37/hg19 with GRCh38."
        )
    if assemblies and ref_assembly and ref_assembly not in assemblies:
        raise ValidationError(
            f"Reference assembly hint '{ref_assembly}' does not match input "
            f"assembly '{next(iter(assemblies))}'."
        )

    resolved = next(iter(assemblies)) if assemblies else (ref_assembly or "unknown")

    # Contig dictionaries among inputs
    base = metas[0]
    for other in metas[1:]:
        ok, problems = contig_dict_compatible(base.contigs, other.contigs)
        if not ok:
            raise ValidationError(
                "Incompatible contig dictionaries between "
                f"'{base.path}' and '{other.path}': " + "; ".join(problems)
            )

    # Sample composition
    sample_sets = [tuple(m.samples) for m in metas]
    if require_identical_samples:
        unique = {frozenset(s) for s in sample_sets}
        if len(unique) > 1:
            detail = ", ".join(f"{m.path}:{list(m.samples)}" for m in metas)
            raise SampleMismatchError(
                "Sample sets differ across caller VCFs; refusing to align by column "
                f"position. Details: {detail}"
            )
        # Also reject differently ordered multi-sample if names match but order differs
        # when lengths > 1 — reordering is OK if names match; we key by name later.
    for m in metas:
        if mode == AnalysisMode.GERMLINE and len(m.samples) == 0:
            raise ValidationError(f"No samples in germline VCF: {m.path}")

    # Mode hints: warn if Mutect2-like headers in germline without somatic mode
    if mode == AnalysisMode.GERMLINE:
        for m in metas:
            lower = os.path.basename(m.path).lower()
            if ".mu." in lower or "mutect" in lower:
                logger.warning(
                    "Input %s looks like Mutect2 while --mode germline; "
                    "somatic semantics will not be applied.",
                    m.path,
                )

    return resolved
