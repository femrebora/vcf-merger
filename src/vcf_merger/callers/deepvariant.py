"""DeepVariant adapter."""

from __future__ import annotations

import re
from typing import Any, Optional

import pysam

from vcf_merger.callers.base import CallerAdapter, _find_version, _header_text


class DeepVariantAdapter(CallerAdapter):
    name = "deepvariant"

    def matches_header(self, header: pysam.VariantHeader) -> tuple[bool, Optional[str]]:
        text = _header_text(header)
        if re.search(r"deepvariant|DeepVariant", text):
            ver = _find_version(
                [re.compile(r"[Dd]eep[Vv]ariant[^0-9]*([0-9]+\.[0-9][^\s\",]*)")],
                text,
            )
            return True, ver
        if "VAF" in header.formats and "PL" in header.formats and "GQ" in header.formats:
            if "source=" in text.lower() and "deepvariant" in text.lower():
                return True, None
        return False, None

    def matches_filename(self, path: str) -> bool:
        p = path.lower()
        return ".dv." in p or "deepvariant" in p

    def caller_metrics(self, record: pysam.VariantRecord, alt_index: int) -> dict[str, Any]:
        return {}
