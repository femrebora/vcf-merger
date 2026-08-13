"""Helpers for computing allele fractions from somatic caller FORMAT fields."""

from __future__ import annotations

from typing import Any, Optional


def _as_float(value: Any) -> Optional[float]:
    if value is None or value == ".":
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        return _as_float(value[0])
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_int(value: Any) -> Optional[int]:
    if value is None or value == ".":
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def af_from_ad(fmt: dict[str, Any], alt_index: int = 0) -> Optional[float]:
    ad = fmt.get("AD")
    if not isinstance(ad, list) or len(ad) < 2:
        return None
    ref = _as_int(ad[0])
    idx = alt_index + 1
    alt = _as_int(ad[idx]) if idx < len(ad) else None
    if ref is None or alt is None:
        return None
    total = ref + alt
    if total <= 0:
        return None
    return alt / total


def strelka_snv_af(fmt: dict[str, Any], alt_base: str) -> Optional[float]:
    """Estimate Strelka SNV allele fraction from tier1 base counts.

    Strelka somatic SNVs typically expose AU,CU,GU,TU as Number=2 (tier1,tier2).
    """
    alt_base = alt_base.upper()
    key = {"A": "AU", "C": "CU", "G": "GU", "T": "TU"}.get(alt_base)
    if not key:
        return None
    counts = {}
    for base, fmt_key in (("A", "AU"), ("C", "CU"), ("G", "GU"), ("T", "TU")):
        raw = fmt.get(fmt_key)
        if raw is None:
            continue
        if isinstance(raw, (list, tuple)):
            counts[base] = _as_int(raw[0]) or 0
        else:
            counts[base] = _as_int(raw) or 0
    if alt_base not in counts:
        return None
    total = sum(counts.values())
    if total <= 0:
        return None
    return counts[alt_base] / total


def pick_sample_fmt(
    samples: dict[str, dict[str, Any]],
    name: Optional[str],
) -> Optional[dict[str, Any]]:
    if not name:
        return None
    return samples.get(name)
