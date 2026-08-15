"""Somatic analysis context (technical only; no AMP/ASCO/CAP tiering)."""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Optional

from vcf_merger.exceptions import ValidationError
from vcf_merger.models import InputVCFMetadata

logger = logging.getLogger(__name__)


@dataclass(slots=True)
class SomaticContext:
    """Tumor / matched-normal / tumor-only roles for somatic mode."""

    tumor_sample: str
    normal_sample: Optional[str] = None
    tumor_only: bool = False

    def describe(self) -> str:
        if self.tumor_only or not self.normal_sample:
            return f"tumor-only (tumor={self.tumor_sample})"
        return f"tumor/normal (tumor={self.tumor_sample}, normal={self.normal_sample})"


def resolve_somatic_context(
    metas: list[InputVCFMetadata],
    *,
    tumor_sample: Optional[str] = None,
    normal_sample: Optional[str] = None,
    tumor_only: bool = False,
) -> SomaticContext:
    """Resolve and validate somatic sample roles."""
    all_samples: list[str] = []
    seen: set[str] = set()
    for m in metas:
        for s in m.samples:
            if s not in seen:
                seen.add(s)
                all_samples.append(s)

    if not all_samples:
        raise ValidationError("Somatic mode requires at least one sample column in input VCFs")

    if tumor_only and normal_sample:
        raise ValidationError("Cannot combine --tumor-only with --normal-sample")

    # Infer tumor when omitted: prefer common names, else first sample
    resolved_tumor = tumor_sample
    if not resolved_tumor:
        for candidate in ("TUMOR", "Tumor", "tumor", "TUMOUR", "Tumour"):
            if candidate in seen:
                resolved_tumor = candidate
                break
        if not resolved_tumor:
            resolved_tumor = all_samples[0]
            logger.warning(
                "Somatic mode: --tumor-sample not set; using '%s' as tumor sample.",
                resolved_tumor,
            )

    if resolved_tumor not in seen:
        raise ValidationError(
            f"--tumor-sample '{resolved_tumor}' not found in inputs; available: {all_samples}"
        )

    resolved_normal = normal_sample
    if tumor_only:
        return SomaticContext(
            tumor_sample=resolved_tumor,
            normal_sample=None,
            tumor_only=True,
        )

    if not resolved_normal:
        # Infer matched normal when exactly two samples and tumor is one of them
        if len(all_samples) == 2:
            resolved_normal = next(s for s in all_samples if s != resolved_tumor)
            logger.info("Inferred matched normal sample '%s'", resolved_normal)
        else:
            for candidate in ("NORMAL", "Normal", "normal", "BLOOD", "Blood"):
                if candidate in seen and candidate != resolved_tumor:
                    resolved_normal = candidate
                    break

    if resolved_normal is not None and resolved_normal not in seen:
        raise ValidationError(
            f"--normal-sample '{resolved_normal}' not found in inputs; available: {all_samples}"
        )

    if resolved_normal is None and len(all_samples) > 1:
        logger.warning(
            "Somatic mode with multiple samples but no --normal-sample; "
            "treating as tumor-only using '%s'. Pass --tumor-only to silence this warning.",
            resolved_tumor,
        )
        return SomaticContext(
            tumor_sample=resolved_tumor,
            normal_sample=None,
            tumor_only=True,
        )

    return SomaticContext(
        tumor_sample=resolved_tumor,
        normal_sample=resolved_normal,
        tumor_only=resolved_normal is None,
    )
