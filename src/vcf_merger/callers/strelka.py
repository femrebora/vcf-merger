"""Strelka / Strelka2 adapter."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

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

    def tumor_af(self, record: pysam.VariantRecord, alt_index: int) -> Optional[float]:
        # Strelka often encodes tier counts in FORMAT AU/CU/GU/TU etc.
        return None
