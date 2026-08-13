"""ID-aware VCF header reconciliation."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any, Optional

import pysam

from vcf_merger.exceptions import HeaderConflictError

logger = logging.getLogger(__name__)

VM_INFO_DEFINITIONS = [
    {
        "id": "VM_CALLERS",
        "number": ".",
        "type": "String",
        "description": "Callers that detected this allele (technical ensemble provenance; not ACMG evidence)",
    },
    {
        "id": "VM_PASS_CALLERS",
        "number": ".",
        "type": "String",
        "description": "Callers with FILTER=PASS for this allele (technical; not clinical significance)",
    },
    {
        "id": "VM_CALLER_COUNT",
        "number": "1",
        "type": "Integer",
        "description": "Number of distinct callers detecting this allele",
    },
    {
        "id": "VM_PASS_CALLER_COUNT",
        "number": "1",
        "type": "Integer",
        "description": "Number of distinct callers with PASS for this allele",
    },
    {
        "id": "VM_GT_CONFLICT",
        "number": "0",
        "type": "Flag",
        "description": "Set when callers disagree on genotype for the representative sample",
    },
    {
        "id": "VM_ORIG",
        "number": "1",
        "type": "String",
        "description": "Original unnormalized representation chrom:pos:ref>alt from a supporting caller",
    },
]


@dataclass
class HeaderField:
    kind: str  # INFO, FORMAT, FILTER, ALT, contig
    id: str
    number: Optional[str] = None
    type: Optional[str] = None
    description: Optional[str] = None
    length: Optional[int] = None
    raw: Optional[str] = None


def _parse_structured(rec: Any) -> Optional[HeaderField]:
    key = getattr(rec, "key", None)
    if key not in ("INFO", "FORMAT", "FILTER", "ALT", "contig"):
        return None
    items = dict(rec.items()) if hasattr(rec, "items") else {}
    fid = items.get("ID") or getattr(rec, "id", None) or items.get("id")
    if key == "contig":
        fid = items.get("ID")
        length = items.get("length")
        return HeaderField(
            kind="contig",
            id=str(fid),
            length=int(length) if length not in (None, ".") else None,
            raw=str(rec),
        )
    return HeaderField(
        kind=key,
        id=str(fid),
        number=str(items["Number"]) if "Number" in items else None,
        type=str(items["Type"]) if "Type" in items else None,
        description=str(items.get("Description", "")).strip('"'),
        raw=str(rec),
    )


def reconcile_headers(
    headers: list[pysam.VariantHeader],
    *,
    samples: list[str],
    namespace_conflicts: bool = True,
    tool_meta_lines: Optional[list[str]] = None,
) -> pysam.VariantHeader:
    """Build a reconciled output header from input headers.

    Conflicting INFO/FORMAT definitions (same ID, different Number/Type) raise
    unless ``namespace_conflicts`` is True, in which case later definitions are
    renamed with a numeric suffix and logged.
    """
    if not headers:
        raise HeaderConflictError("No headers to reconcile")

    out = pysam.VariantHeader()
    out.add_line("##fileformat=VCFv4.2")

    seen_generic: set[str] = set()
    fields: dict[tuple[str, str], HeaderField] = {}

    for header in headers:
        for rec in header.records:
            key = getattr(rec, "key", "")
            if key == "fileformat":
                continue
            parsed = _parse_structured(rec)
            if parsed is None:
                line = str(rec).rstrip("\n")
                if not line.startswith("##"):
                    line = "##" + line
                if line not in seen_generic:
                    try:
                        out.add_line(line)
                        seen_generic.add(line)
                    except Exception:  # noqa: BLE001
                        logger.debug("Skipping unparsable header line: %s", line)
                continue

            map_key = (parsed.kind, parsed.id)
            if map_key not in fields:
                fields[map_key] = parsed
                continue

            existing = fields[map_key]
            if parsed.kind in ("INFO", "FORMAT"):
                if (existing.number, existing.type) != (parsed.number, parsed.type):
                    if not namespace_conflicts:
                        raise HeaderConflictError(
                            f"Conflicting {parsed.kind} definition for ID={parsed.id}: "
                            f"Number/Type {existing.number}/{existing.type} vs "
                            f"{parsed.number}/{parsed.type}"
                        )
                    # Namespace the conflicting definition
                    suffix = 2
                    new_id = f"{parsed.id}_VM{suffix}"
                    while (parsed.kind, new_id) in fields:
                        suffix += 1
                        new_id = f"{parsed.id}_VM{suffix}"
                    logger.warning(
                        "Namespacing conflicting %s ID %s -> %s",
                        parsed.kind,
                        parsed.id,
                        new_id,
                    )
                    parsed.id = new_id
                    fields[(parsed.kind, new_id)] = parsed
            elif parsed.kind == "contig":
                if (
                    existing.length is not None
                    and parsed.length is not None
                    and existing.length != parsed.length
                ):
                    raise HeaderConflictError(
                        f"Conflicting contig length for {parsed.id}: "
                        f"{existing.length} vs {parsed.length}"
                    )

    # Emit contigs in first-seen order from first header preferentially
    for header in headers:
        for name in header.contigs:
            key = ("contig", name)
            if key in fields:
                f = fields.pop(key)
                length_part = f",length={f.length}" if f.length else ""
                out.add_line(f"##contig=<ID={f.id}{length_part}>")

    for (kind, fid), f in list(fields.items()):
        if kind == "contig":
            length_part = f",length={f.length}" if f.length else ""
            out.add_line(f"##contig=<ID={f.id}{length_part}>")
            fields.pop((kind, fid), None)

    for (kind, fid), f in fields.items():
        if kind == "INFO":
            out.add_line(
                f'##INFO=<ID={f.id},Number={f.number or "."},Type={f.type or "String"},'
                f'Description="{f.description or ""}">'
            )
        elif kind == "FORMAT":
            out.add_line(
                f'##FORMAT=<ID={f.id},Number={f.number or "."},Type={f.type or "String"},'
                f'Description="{f.description or ""}">'
            )
        elif kind == "FILTER":
            out.add_line(
                f'##FILTER=<ID={f.id},Description="{f.description or ""}">'
            )
        elif kind == "ALT":
            out.add_line(
                f'##ALT=<ID={f.id},Description="{f.description or ""}">'
            )

    # Ensure PASS filter exists
    if "PASS" not in out.filters:
        out.add_line('##FILTER=<ID=PASS,Description="All filters passed">')

    # Add VM_* fields if not conflicting
    existing_info = set(out.info.keys())
    for spec in VM_INFO_DEFINITIONS:
        if spec["id"] in existing_info:
            raise HeaderConflictError(
                f"Output INFO ID {spec['id']} already defined in inputs; "
                f"refusing to overwrite. Rename or remove the conflicting definition."
            )
        out.add_line(
            f"##INFO=<ID={spec['id']},Number={spec['number']},Type={spec['type']},"
            f"Description=\"{spec['description']}\">"
        )

    # Minimal FORMAT for representative genotype
    if "GT" not in out.formats:
        out.add_line('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    if "DP" not in out.formats:
        out.add_line('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    if "GQ" not in out.formats:
        out.add_line(
            '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">'
        )

    if tool_meta_lines:
        for line in tool_meta_lines:
            if not line.startswith("##"):
                line = "##" + line
            out.add_line(line)

    for sample in samples:
        out.add_sample(sample)

    return out
