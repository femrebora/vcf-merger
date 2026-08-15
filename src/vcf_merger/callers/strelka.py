"""Strelka / Strelka2 adapter."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

from vcf_merger.af_utils import af_from_ad, strelka_snv_af
from vcf_merger.callers.base import CallerAdapter, _as_float, _find_version, _header_text


class StrelkaAdapter(CallerAdapter):
    name = "strelka"

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        text = _header_text(header)
        if re.search(r"strelka", text, re.IGNORECASE):
            ver = _find_version(
                [re.compile(r"[Ss]trelka(?:2)?[^0-9]*([0-9]+\.[0-9][^\s\",]*)")],
                text,
            )
            return True, ver
        # somatic SNV fields
        if "SomaticEVS" in header.info or "QSS_NT" in header.info or "QSI_NT" in header.info:
            return True, None
        return False, None

    def matches_filename(self, path: str) -> bool:
        p = path.lower()
        return ".st." in p or "strelka" in p

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        metrics: dict[str, Any] = {}
        for key in ("SomaticEVS", "QSS_NT", "QSI_NT", "TQSS", "NT", "SGT"):
            if key in record.info:
                metrics[key] = record.info[key]
        return metrics

    def _af_for_sample(
        self,
        record: pysam.VariantRecord,
        alt_index: int,
        sample_fmts: Optional[dict[str, dict[str, Any]]],
        sample: Optional[str],
    ) -> Optional[float]:
        if not sample_fmts or not sample or sample not in sample_fmts:
            return None
        fmt = sample_fmts[sample]
        alt = (record.alts or (".",))[alt_index]
        af = af_from_ad(fmt, alt_index)
        if af is not None:
            return af
        if isinstance(alt, str) and len(alt) == 1:
            return strelka_snv_af(fmt, alt)
        return None

    def tumor_af(
        self,
        record: pysam.VariantRecord,
        alt_index: int,
        *,
        sample_fmts: Optional[dict[str, dict[str, Any]]] = None,
        tumor_sample: Optional[str] = None,
    ) -> Optional[float]:
        return self._af_for_sample(record, alt_index, sample_fmts, tumor_sample)

    def normal_af(
        self,
        record: pysam.VariantRecord,
        alt_index: int,
        *,
        sample_fmts: Optional[dict[str, dict[str, Any]]] = None,
        normal_sample: Optional[str] = None,
    ) -> Optional[float]:
        return self._af_for_sample(record, alt_index, sample_fmts, normal_sample)
