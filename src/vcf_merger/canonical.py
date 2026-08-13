"""Canonical small-variant identity helpers."""

from __future__ import annotations

from typing import Optional

from vcf_merger.models import CanonicalVariant, VrsIdentifier


def make_canonical(
    assembly: str,
    contig: str,
    pos: int,
    ref: str,
    alt: str,
    *,
    vrs: Optional[VrsIdentifier] = None,
) -> CanonicalVariant:
    """Build a canonical identity from normalized alleles.

    Does not use raw unnormalized CHROM/POS/REF/ALT alone as identity.
    """
    assembly_n = (assembly or "unknown").strip()
    contig_n = contig.strip()
    ref_n = ref.upper()
    alt_n = alt.upper()
    vrs_id = None
    if vrs is not None:
        vrs_id = vrs.identify(assembly_n, contig_n, pos, ref_n, alt_n)
    return CanonicalVariant(
        assembly=assembly_n,
        contig=contig_n,
        pos=int(pos),
        ref=ref_n,
        alt=alt_n,
        vrs_id=vrs_id,
    )


def is_symbolic_allele(allele: str) -> bool:
    """Return True for symbolic / breakend-style alleles."""
    if not allele:
        return False
    if allele.startswith("<") and allele.endswith(">"):
        return True
    if "[" in allele or "]" in allele:
        return True  # BND
    if allele == "*":
        return True
    return False


def is_supported_small_variant(ref: str, alt: str) -> tuple[bool, Optional[str]]:
    """First-class support: SNV / small indel / MNV after normalization.

    Structural alleles are rejected for the small-variant harmonizer.
    """
    if is_symbolic_allele(alt) or is_symbolic_allele(ref):
        return False, f"symbolic_or_sv_allele:{alt}"
    if not ref or not alt:
        return False, "empty_allele"
    # Allow ACGT and IUPAC for small variants; reject clearly structural tokens
    allowed = set("ACGTNacgtn")
    if not set(ref) <= allowed or not set(alt) <= allowed:
        return False, f"non_acgt_allele:{ref}>{alt}"
    return True, None
