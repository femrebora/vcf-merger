"""GATK HaplotypeCaller adapter."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

from vcf_merger.callers.base import CallerAdapter, _find_version, _header_text


class HaplotypeCallerAdapter(CallerAdapter):
    name = "haplotypecaller"

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        text = _header_text(header)
        if re.search(r"HaplotypeCaller", text, re.IGNORECASE):
            ver = _find_version(
                [
                    re.compile(r"HaplotypeCaller.*?Version[=/ ]([0-9][^,\s\"]*)", re.I),
                    re.compile(r"GATK.*?([0-9]+\.[0-9][^\s\",]*)", re.I),
                ],
                text,
            )
            return True, ver
        # Common HC FORMAT fields without Mutect2 markers
        if "HaplotypeCaller" not in text and "Mutect2" not in text:
            if "F1R2" in header.formats and "F2R1" in header.formats and "AD" in header.formats:
                # weak signal — defer to filename usually
                return False, None
        return False, None

    def matches_filename(self, path: str) -> bool:
        p = path.lower()
        return ".hc." in p or "haplotypecaller" in p

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        metrics: dict[str, Any] = {}
        for key in ("MQ", "MQRankSum", "ReadPosRankSum", "BaseQRankSum", "ClippingRankSum", "DP"):
            if key in record.info:
                metrics[key] = record.info[key]
        return metrics
