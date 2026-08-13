"""FreeBayes adapter."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

from vcf_merger.callers.base import CallerAdapter, _find_version, _header_text


class FreeBayesAdapter(CallerAdapter):
    name = "freebayes"

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        text = _header_text(header)
        if re.search(r"freebayes", text, re.IGNORECASE):
            ver = _find_version(
                [re.compile(r"freeBayes\s+([0-9][^,\s\"]*)", re.I), re.compile(r"freebayes\s+([0-9][^,\s\"]*)", re.I)],
                text,
            )
            return True, ver
        # characteristic INFO fields
        if "AO" in header.info and "RO" in header.info and "AB" in header.info:
            return True, None
        return False, None

    def matches_filename(self, path: str) -> bool:
        p = path.lower()
        return ".fb." in p or "freebayes" in p

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        metrics: dict[str, Any] = {}
        for key in ("AO", "RO", "AB", "SAF", "SAR", "SRF", "SRR", "DP"):
            if key in record.info:
                val = record.info[key]
                if isinstance(val, tuple) and len(val) > alt_index:
                    metrics[key] = val[alt_index]
                else:
                    metrics[key] = val
        return metrics
