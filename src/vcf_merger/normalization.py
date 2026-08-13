"""Reference-aware normalization via bcftools norm."""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Optional

from vcf_merger.exceptions import ExternalToolError, NormalizationError, ReferenceError
from vcf_merger.validation import validate_reference_fasta

logger = logging.getLogger(__name__)


def require_bcftools() -> str:
    path = shutil.which("bcftools")
    if not path:
        raise ExternalToolError(
            "bcftools not found on PATH. Install bcftools to run normalization "
            "(https://samtools.github.io/bcftools/)."
        )
    return path


def normalize_vcf(
    input_path: str | Path,
    output_path: str | Path,
    *,
    reference: str | Path,
    multiallelic: str = "-any",
    check_ref: str = "w",
    threads: int = 1,
) -> str:
    """Normalize a VCF with bcftools norm (left-align, trim, optional split).

    Parameters
    ----------
    multiallelic :
        Passed to ``bcftools norm -m``. Default ``-any`` splits multiallelics.
    check_ref :
        Passed to ``bcftools norm -c`` (e.g. ``w`` warn, ``e`` error, ``x`` exclude).

    Returns
    -------
    str
        Path to the normalized output.
    """
    bcftools = require_bcftools()
    ref = validate_reference_fasta(reference)
    in_path = str(input_path)
    out_path = str(output_path)

    if not os.path.exists(in_path):
        raise NormalizationError(f"Input VCF not found: {in_path}")

    out_dir = os.path.dirname(os.path.abspath(out_path)) or "."
    os.makedirs(out_dir, exist_ok=True)

    cmd = [
        bcftools,
        "norm",
        "-f",
        ref,
        "-m",
        multiallelic,
        "-c",
        check_ref,
        "--threads",
        str(threads),
        "-O",
        "z" if out_path.endswith(".gz") else "v",
        "-o",
        out_path,
        in_path,
    ]
    logger.info("Running: %s", " ".join(cmd))
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    except OSError as exc:
        raise ExternalToolError(f"Failed to execute bcftools: {exc}") from exc

    if proc.returncode != 0:
        raise NormalizationError(
            f"bcftools norm failed (exit {proc.returncode}): {proc.stderr.strip()}"
        )
    if proc.stderr.strip():
        logger.debug("bcftools norm stderr: %s", proc.stderr.strip())

    # Index gz outputs when tabix available
    if out_path.endswith(".gz"):
        tabix = shutil.which("tabix")
        if tabix:
            subprocess.run([tabix, "-f", "-p", "vcf", out_path], check=False)

    return out_path


def normalize_to_temp(
    input_path: str | Path,
    *,
    reference: str | Path,
    multiallelic: str = "-any",
    check_ref: str = "w",
) -> str:
    """Normalize into a temporary .vcf.gz; caller must delete when done."""
    fd, tmp = tempfile.mkstemp(prefix="vcf_merger_norm_", suffix=".vcf.gz")
    os.close(fd)
    try:
        return normalize_vcf(
            input_path,
            tmp,
            reference=reference,
            multiallelic=multiallelic,
            check_ref=check_ref,
        )
    except Exception:
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
