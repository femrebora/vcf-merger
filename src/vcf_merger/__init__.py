"""vcf-merger: technical VCF harmonization and ensemble evidence aggregation.

This package does not perform ACMG/AMP pathogenicity classification or
clinical interpretation. Caller concordance is a technical metric only.
"""

from __future__ import annotations

__version__ = "0.2.0"

from vcf_merger.harmonizer import harmonize_vcfs
from vcf_merger.inspect import inspect_vcf
from vcf_merger.normalization import normalize_vcf

__all__ = [
    "__version__",
    "harmonize_vcfs",
    "inspect_vcf",
    "normalize_vcf",
]
