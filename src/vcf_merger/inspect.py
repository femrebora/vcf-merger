"""Inspect VCF/BCF files and report metadata."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any, Optional

from vcf_merger.callers import detect_caller
from vcf_merger.models import DetectionSource, InputVCFMetadata
from vcf_merger.validation import collect_metadata


def inspect_vcf(
    path: str | Path,
    *,
    caller: Optional[str] = None,
) -> InputVCFMetadata:
    """Return metadata for one VCF, including layered caller detection."""
    meta = collect_metadata(path)
    detected = detect_caller(path, header_meta=meta, explicit=caller)
    meta.caller = detected.name
    meta.caller_version = detected.version
    meta.caller_detection_source = detected.source
    return meta


def metadata_to_dict(meta: InputVCFMetadata) -> dict[str, Any]:
    return {
        "path": meta.path,
        "file_kind": meta.file_kind.value,
        "fileformat": meta.fileformat,
        "assembly": meta.assembly,
        "reference": meta.reference,
        "samples": meta.samples,
        "contigs": [{"name": c.name, "length": c.length} for c in meta.contigs],
        "caller": meta.caller,
        "caller_version": meta.caller_version,
        "caller_detection_source": (
            meta.caller_detection_source.value
            if isinstance(meta.caller_detection_source, DetectionSource)
            else str(meta.caller_detection_source)
        ),
        "is_gvcf": meta.is_gvcf,
        "gvcf_reasons": meta.gvcf_reasons,
        "info_ids": meta.header_info_ids,
        "format_ids": meta.header_format_ids,
        "filter_ids": meta.header_filter_ids,
    }


def inspect_vcf_json(path: str | Path, *, caller: Optional[str] = None) -> str:
    return json.dumps(metadata_to_dict(inspect_vcf(path, caller=caller)), indent=2)
