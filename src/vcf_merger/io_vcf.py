"""VCF I/O helpers built on pysam.VariantFile."""

from __future__ import annotations

import logging
import os
from collections.abc import Iterator
from pathlib import Path
from typing import Any, Optional

import pysam

from vcf_merger.exceptions import ValidationError
from vcf_merger.models import ContigInfo, FileKind

logger = logging.getLogger(__name__)


def detect_file_kind(path: str | Path) -> FileKind:
    name = os.path.basename(str(path)).lower()
    if name.endswith(".g.vcf.gz") or name.endswith(".gvcf.gz"):
        return FileKind.GVCF_GZ
    if name.endswith(".g.vcf") or name.endswith(".gvcf"):
        return FileKind.GVCF
    if name.endswith(".bcf"):
        return FileKind.BCF
    if name.endswith(".vcf.gz"):
        return FileKind.VCF_GZ
    if name.endswith(".vcf"):
        return FileKind.VCF
    return FileKind.UNKNOWN


def open_variant_file(path: str | Path, mode: str = "r") -> pysam.VariantFile:
    """Open a VCF/BCF with pysam."""
    path_s = str(path)
    if not os.path.exists(path_s) and "w" not in mode:
        raise ValidationError(f"VCF not found: {path_s}")
    try:
        return pysam.VariantFile(path_s, mode)
    except Exception as exc:  # noqa: BLE001 - surface as domain error
        raise ValidationError(f"Failed to open VCF '{path_s}': {exc}") from exc


def iter_records(path: str | Path) -> Iterator[pysam.VariantRecord]:
    with open_variant_file(path) as vf:
        yield from vf


def contig_list_from_header(header: pysam.VariantHeader) -> list[ContigInfo]:
    contigs: list[ContigInfo] = []
    for name in header.contigs:
        ctg = header.contigs[name]
        length = getattr(ctg, "length", None)
        contigs.append(ContigInfo(name=name, length=length if length and length > 0 else None))
    return contigs


def header_assembly_hints(header: pysam.VariantHeader) -> tuple[Optional[str], Optional[str]]:
    """Best-effort assembly and reference metadata from header."""
    reference = None
    assembly = None

    # ##reference=...
    for rec in header.records:
        key = rec.key.lower() if hasattr(rec, "key") else ""
        if key == "reference":
            reference = rec.value if hasattr(rec, "value") else str(rec)
        if key == "assembly":
            assembly = rec.value if hasattr(rec, "value") else str(rec)
        # contig assembly attribute
        if key == "contig" and assembly is None:
            items = dict(rec.items()) if hasattr(rec, "items") else {}
            if "assembly" in items:
                assembly = items["assembly"]

    text = " ".join(
        filter(
            None,
            [reference or "", assembly or "", str(header)],
        )
    ).lower()

    if assembly is None:
        if "grch38" in text or "hg38" in text or "hs38" in text:
            assembly = "GRCh38"
        elif "grch37" in text or "hg19" in text or "hs37" in text:
            assembly = "GRCh37"

    return assembly, reference


def record_to_info_dict(record: pysam.VariantRecord) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for key, value in record.info.items():
        out[key] = _jsonable(value)
    return out


def sample_format_dict(record: pysam.VariantRecord, sample: str) -> dict[str, Any]:
    data = record.samples[sample]
    out: dict[str, Any] = {}
    for key in data.keys():
        try:
            out[key] = _jsonable(data[key])
        except Exception:  # noqa: BLE001
            out[key] = None
    return out


def genotype_string(record: pysam.VariantRecord, sample: str) -> tuple[Optional[str], bool, Optional[int]]:
    """Return (GT string, phased, ploidy) preserving source ploidy."""
    data = record.samples[sample]
    alleles = data.get("GT")
    if alleles is None:
        return None, False, None
    # pysam may return tuple of ints/None
    if not isinstance(alleles, (list, tuple)):
        return None, bool(data.phased), None
    phased = bool(data.phased)
    sep = "|" if phased else "/"
    parts: list[str] = []
    for a in alleles:
        if a is None:
            parts.append(".")
        else:
            parts.append(str(a))
    if not parts:
        return None, phased, 0
    return sep.join(parts), phased, len(parts)


def _jsonable(value: Any) -> Any:
    if value is None:
        return None
    if isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="replace")
    if isinstance(value, tuple):
        return [_jsonable(v) for v in value]
    if isinstance(value, list):
        return [_jsonable(v) for v in value]
    return str(value)


def index_variant_file(path: str | Path) -> str:
    """Create tabix or CSI index for a bgzipped VCF. Returns index path."""
    path_s = str(path)
    pysam.tabix_index(path_s, preset="vcf", force=True)
    tbi = path_s + ".tbi"
    csi = path_s + ".csi"
    if os.path.exists(tbi):
        return tbi
    if os.path.exists(csi):
        return csi
    raise ValidationError(f"Failed to create index for {path_s}")
