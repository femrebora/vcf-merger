"""Caller adapter base types and layered detection."""

from __future__ import annotations

import os
import re
from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional

import pysam

from vcf_merger.io_vcf import genotype_string, record_to_info_dict, sample_format_dict
from vcf_merger.models import (
    CallerEvidence,
    DetectionSource,
    InputVCFMetadata,
    SampleEvidence,
)


@dataclass(slots=True)
class DetectedCaller:
    name: str
    version: Optional[str]
    source: DetectionSource
    confidence: str  # high | medium | low


class CallerAdapter(ABC):
    """Extract technical evidence without forcing identical native fields."""

    name: str

    @abstractmethod
    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        """Return (matched, version)."""

    def matches_filename(self, path: str) -> bool:
        return False

    def extract_evidence(
        self,
        record: pysam.VariantRecord,
        *,
        alt_index: int,
        source_file: str,
        caller_version: Optional[str],
        detection_source: DetectionSource,
        tumor_sample: Optional[str] = None,
        normal_sample: Optional[str] = None,
    ) -> CallerEvidence:
        alt = (record.alts or (".",))[alt_index]
        filt = list(record.filter.keys()) if record.filter else []
        if not filt:
            filter_status = "."
            passed = True  # missing FILTER treated cautiously as unfiltered
        else:
            filter_status = ";".join(filt)
            passed = filt == ["PASS"] or (len(filt) == 1 and filt[0] == "PASS")

        samples: list[SampleEvidence] = []
        sample_fmts: dict[str, dict[str, Any]] = {}
        for sample in record.samples:
            gt, phased, ploidy = genotype_string(record, sample)
            fmt = sample_format_dict(record, sample)
            sample_fmts[sample] = fmt
            depth = _as_int(fmt.get("DP"))
            ref_depth, alt_depth, af = self._depths_and_af(fmt, alt_index)
            gq = _as_float(fmt.get("GQ"))
            samples.append(
                SampleEvidence(
                    sample=sample,
                    genotype=gt,
                    phased=phased,
                    ploidy=ploidy,
                    depth=depth,
                    ref_depth=ref_depth,
                    alt_depth=alt_depth,
                    allele_fraction=af,
                    genotype_quality=gq,
                    raw_format=fmt,
                )
            )

        return CallerEvidence(
            caller=self.name,
            caller_version=caller_version,
            source_file=source_file,
            filter_status=filter_status,
            passed=passed,
            qual=_as_float(record.qual),
            samples=samples,
            caller_specific_metrics=self.caller_metrics(record, alt_index),
            raw_info=record_to_info_dict(record),
            original_chrom=str(record.contig),
            original_pos=int(record.pos),
            original_ref=str(record.ref),
            original_alt=str(alt),
            original_record_identifier=f"{record.contig}:{record.pos}:{record.ref}>{alt}",
            detection_source=detection_source,
            tumor_allele_fraction=self.tumor_af(
                record, alt_index, sample_fmts=sample_fmts, tumor_sample=tumor_sample
            ),
            normal_allele_fraction=self.normal_af(
                record, alt_index, sample_fmts=sample_fmts, normal_sample=normal_sample
            ),
        )

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        return {}

    def tumor_af(
        self,
        record: pysam.VariantRecord,
        alt_index: int,
        *,
        sample_fmts: Optional[dict[str, dict[str, Any]]] = None,
        tumor_sample: Optional[str] = None,
    ) -> Optional[float]:
        return None

    def normal_af(
        self,
        record: pysam.VariantRecord,
        alt_index: int,
        *,
        sample_fmts: Optional[dict[str, dict[str, Any]]] = None,
        normal_sample: Optional[str] = None,
    ) -> Optional[float]:
        return None

    def _depths_and_af(
        self, fmt: dict[str, Any], alt_index: int
    ) -> tuple[Optional[int], Optional[int], Optional[float]]:
        ad = fmt.get("AD")
        ref_depth = alt_depth = af = None
        if isinstance(ad, list) and len(ad) >= 2:
            ref_depth = _as_int(ad[0])
            idx = alt_index + 1
            if idx < len(ad):
                alt_depth = _as_int(ad[idx])
            if ref_depth is not None and alt_depth is not None and (ref_depth + alt_depth) > 0:
                af = alt_depth / (ref_depth + alt_depth)
        if af is None:
            for key in ("AF", "VAF", "FREQ"):
                if key in fmt:
                    af = _as_float(fmt[key] if not isinstance(fmt[key], list) else fmt[key][0])
                    break
        return ref_depth, alt_depth, af


def _as_int(value: Any) -> Optional[int]:
    if value is None or value == ".":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _as_float(value: Any) -> Optional[float]:
    if value is None or value == ".":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _header_text(header: pysam.VariantHeader) -> str:
    return str(header)


def _find_version(patterns: list[re.Pattern[str]], text: str) -> Optional[str]:
    for pat in patterns:
        m = pat.search(text)
        if m:
            return m.group(1)
    return None


# Filename fallback tags (not primary)
_FILENAME_TAGS: list[tuple[str, str]] = [
    (".mu.", "mutect2"),
    (".fb.", "freebayes"),
    (".hc.", "haplotypecaller"),
    (".st.", "strelka"),
    (".dv.", "deepvariant"),
    ("mutect", "mutect2"),
    ("freebayes", "freebayes"),
    ("haplotypecaller", "haplotypecaller"),
    ("deepvariant", "deepvariant"),
    ("strelka", "strelka"),
]


def filename_caller_guess(path: str | Path) -> Optional[str]:
    base = os.path.basename(str(path)).lower()
    for tag, name in _FILENAME_TAGS:
        if tag in base:
            return name
    return None
