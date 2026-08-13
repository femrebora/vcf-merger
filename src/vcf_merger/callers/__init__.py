"""Caller registry and layered detection."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import pysam

from vcf_merger.callers.base import CallerAdapter, DetectedCaller, filename_caller_guess
from vcf_merger.callers.deepvariant import DeepVariantAdapter
from vcf_merger.callers.freebayes import FreeBayesAdapter
from vcf_merger.callers.haplotypecaller import HaplotypeCallerAdapter
from vcf_merger.callers.mutect2 import Mutect2Adapter
from vcf_merger.callers.strelka import StrelkaAdapter
from vcf_merger.io_vcf import open_variant_file
from vcf_merger.models import DetectionSource, InputVCFMetadata

_ADAPTERS: list[CallerAdapter] = [
    Mutect2Adapter(),
    FreeBayesAdapter(),
    HaplotypeCallerAdapter(),
    StrelkaAdapter(),
    DeepVariantAdapter(),
]

_BY_NAME = {a.name: a for a in _ADAPTERS}

# Short aliases used in legacy CALLERS tags / CLI
ALIASES = {
    "mu": "mutect2",
    "mutect2": "mutect2",
    "mutect": "mutect2",
    "fb": "freebayes",
    "freebayes": "freebayes",
    "hc": "haplotypecaller",
    "haplotypecaller": "haplotypecaller",
    "st": "strelka",
    "strelka": "strelka",
    "strelka2": "strelka",
    "dv": "deepvariant",
    "deepvariant": "deepvariant",
    "unknown": "unknown",
}


def normalize_caller_name(name: str) -> str:
    return ALIASES.get(name.strip().lower(), name.strip().lower())


def get_adapter(name: str) -> Optional[CallerAdapter]:
    return _BY_NAME.get(normalize_caller_name(name))


def all_adapters() -> list[CallerAdapter]:
    return list(_ADAPTERS)


def detect_caller(
    path: str | Path,
    *,
    header_meta: Optional[InputVCFMetadata] = None,
    explicit: Optional[str] = None,
) -> DetectedCaller:
    """Layered detection: CLI > header metadata > filename heuristic."""
    if explicit:
        name = normalize_caller_name(explicit)
        return DetectedCaller(name=name, version=None, source=DetectionSource.CLI, confidence="high")

    path_s = str(path)
    with open_variant_file(path_s) as vf:
        header = vf.header
        for adapter in _ADAPTERS:
            matched, version = adapter.matches_header(header)
            if matched:
                return DetectedCaller(
                    name=adapter.name,
                    version=version,
                    source=DetectionSource.HEADER,
                    confidence="high" if version else "medium",
                )

    guess = filename_caller_guess(path_s)
    if guess:
        return DetectedCaller(
            name=guess,
            version=None,
            source=DetectionSource.FILENAME,
            confidence="low",
        )
    return DetectedCaller(
        name="unknown",
        version=None,
        source=DetectionSource.UNKNOWN,
        confidence="low",
    )


def short_label(caller_name: str) -> str:
    """Compact label for INFO fields."""
    mapping = {
        "mutect2": "MU",
        "freebayes": "FB",
        "haplotypecaller": "HC",
        "strelka": "ST",
        "deepvariant": "DV",
        "unknown": "UNK",
    }
    return mapping.get(normalize_caller_name(caller_name), caller_name.upper()[:4])
