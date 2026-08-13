"""GATK Mutect2 adapter (somatic-oriented evidence)."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

from vcf_merger.callers.base import CallerAdapter, _as_float, _find_version, _header_text


class Mutect2Adapter(CallerAdapter):
    name = "mutect2"

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        text = _header_text(header)
        if re.search(r"Mutect2|MUTECT", text):
            ver = _find_version(
                [
                    re.compile(r"Mutect2.*?Version[=/ ]([0-9][^,\s\"]*)", re.I),
                    re.compile(r"GATK.*?([0-9]+\.[0-9][^\s\",]*)", re.I),
                ],
                text,
            )
            return True, ver
        if "TLOD" in header.info or "NLOD" in header.info:
            return True, None
        return False, None

    def matches_filename(self, path: str) -> bool:
        p = path.lower()
        return ".mu." in p or "mutect" in p

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        metrics: dict[str, Any] = {}
        for key in ("TLOD", "NLOD", "GERMQ", "MPOS", "MMQ", "MBQ", "POPAF"):
            if key in record.info:
                val = record.info[key]
                if isinstance(val, tuple) and len(val) > alt_index:
                    metrics[key] = val[alt_index]
                else:
                    metrics[key] = val
        return metrics

    def tumor_af(self, record: pysam.VariantRecord, alt_index: int) -> Optional[float]:
        if "AF" in record.info:
            val = record.info["AF"]
            if isinstance(val, tuple) and len(val) > alt_index:
                return _as_float(val[alt_index])
            return _as_float(val)
        return None
