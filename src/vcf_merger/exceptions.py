"""Domain-specific exceptions for vcf-merger."""

from __future__ import annotations


class VcfMergerError(Exception):
    """Base error for vcf-merger."""


class ValidationError(VcfMergerError):
    """Input VCFs or reference fail validation."""


class GVcfRejectedError(ValidationError):
    """gVCF / reference-confidence input rejected for ensemble merging."""


class ReferenceError(ValidationError):
    """Reference FASTA missing, unusable, or incompatible."""


class HeaderConflictError(ValidationError):
    """Conflicting structured VCF header definitions."""


class SampleMismatchError(ValidationError):
    """Incompatible sample sets across inputs."""


class UnsupportedVariantError(VcfMergerError):
    """Variant class not supported by the small-variant harmonizer."""


class NormalizationError(VcfMergerError):
    """bcftools norm or normalization pipeline failed."""


class CallerDetectionError(VcfMergerError):
    """Unable to determine or configure a variant caller."""


class ExternalToolError(VcfMergerError):
    """Required external tool missing or failed."""
