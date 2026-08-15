"""Machine-readable provenance records."""

from __future__ import annotations

import hashlib
import json
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

from vcf_merger import __version__
from vcf_merger.models import InputVCFMetadata, ProvenanceRecord


def file_checksum(path: str | Path, *, max_bytes: Optional[int] = None) -> str:
    """SHA256 of file contents (full file unless max_bytes set)."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        if max_bytes is None:
            for chunk in iter(lambda: fh.read(1024 * 1024), b""):
                h.update(chunk)
        else:
            h.update(fh.read(max_bytes))
    return h.hexdigest()


def build_provenance(
    *,
    command_line: list[str],
    analysis_mode: str,
    ensemble_strategy: str,
    reference_path: Optional[str],
    normalization: dict[str, Any],
    input_metas: list[InputVCFMetadata],
    samples: list[str],
    notes: Optional[list[str]] = None,
) -> ProvenanceRecord:
    inputs: list[dict[str, Any]] = []
    for m in input_metas:
        entry: dict[str, Any] = {
            "path": m.path,
            "checksum_sha256": file_checksum(m.path) if os.path.exists(m.path) else None,
            "caller": m.caller,
            "caller_version": m.caller_version,
            "caller_detection_source": getattr(m.caller_detection_source, "value", m.caller_detection_source),
            "assembly": m.assembly,
            "samples": m.samples,
            "file_kind": m.file_kind.value if hasattr(m.file_kind, "value") else m.file_kind,
        }
        inputs.append(entry)

    ref_checksum = None
    if reference_path and os.path.exists(reference_path):
        # Reference FASTA can be huge; hash first 64MB + size for practicality note
        size = os.path.getsize(reference_path)
        if size <= 64 * 1024 * 1024:
            ref_checksum = file_checksum(reference_path)
        else:
            ref_checksum = file_checksum(reference_path, max_bytes=64 * 1024 * 1024) + f":partial64MiB:size={size}"

    return ProvenanceRecord(
        tool_name="vcf-merger",
        tool_version=__version__,
        timestamp=datetime.now(timezone.utc).isoformat(),
        command_line=command_line,
        analysis_mode=analysis_mode,
        ensemble_strategy=ensemble_strategy,
        reference_path=reference_path,
        reference_checksum=ref_checksum,
        normalization=normalization,
        inputs=inputs,
        samples=samples,
        notes=notes or [
            "Technical VCF harmonization only; not ACMG/AMP classification.",
            "Caller concordance is not clinical evidence strength.",
        ],
    )


def provenance_to_dict(prov: ProvenanceRecord) -> dict[str, Any]:
    return {
        "tool_name": prov.tool_name,
        "tool_version": prov.tool_version,
        "timestamp": prov.timestamp,
        "command_line": prov.command_line,
        "analysis_mode": prov.analysis_mode,
        "ensemble_strategy": prov.ensemble_strategy,
        "reference_path": prov.reference_path,
        "reference_checksum": prov.reference_checksum,
        "normalization": prov.normalization,
        "inputs": prov.inputs,
        "samples": prov.samples,
        "notes": prov.notes,
    }


def write_provenance(path: str | Path, prov: ProvenanceRecord) -> str:
    path_s = str(path)
    with open(path_s, "w", encoding="utf-8") as fh:
        json.dump(provenance_to_dict(prov), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return path_s
