# vcf_utils.py
#
# Legacy compatibility shim. New code should import from vcf_merger.

from __future__ import annotations

import warnings
from typing import Optional

from vcf_merger.harmonizer import harmonize_vcfs

warnings.warn(
    "vcf_utils is deprecated; import from vcf_merger instead. See docs/migration.md.",
    DeprecationWarning,
    stacklevel=2,
)

# Historical constants retained only for import compatibility.
MERGE_KEY = ["#CHROM", "POS", "REF", "ALT"]
PRIORITY_ORDER = [".MU.vcf", ".FB.vcf", ".HC.vcf", ".ST.vcf", ".DV.vcf"]


def sort_vcf_files(vcf_files: list[str]) -> list[str]:
    """Deprecated no-op sort retained for import compatibility."""
    return list(vcf_files)


def merge_vcfs(
    vcf_files: list[str],
    output_file: str,
    *,
    reference: Optional[str] = None,
    normalize: bool = True,
) -> None:
    """Deprecated wrapper around harmonize_vcfs (union strategy)."""
    if normalize and not reference:
        raise ValueError(
            "merge_vcfs now requires reference=... for normalization, "
            "or normalize=False (not recommended)."
        )
    harmonize_vcfs(
        vcf_files,
        output_file,
        reference=reference,
        normalize=normalize,
        mode="germline",
        strategy="union",
        command_line=["vcf_utils.merge_vcfs"],
    )
