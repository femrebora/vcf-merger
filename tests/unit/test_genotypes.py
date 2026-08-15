"""Additional genotype / allele edge cases."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from vcf_merger.harmonizer import harmonize_vcfs
from vcf_merger.io_vcf import genotype_string, open_variant_file


def test_phased_genotype_preserved_in_evidence(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "phased.vcf.gz"
    result = harmonize_vcfs([fixtures_dir / "phased.HC.vcf"], out, reference=mini_ref)
    with open(result["paths"]["evidence"], encoding="utf-8") as fh:
        row = json.loads(fh.readline())
    assert row["evidence"][0]["samples"][0]["genotype"] == "0|1"
    assert row["evidence"][0]["samples"][0]["phased"] is True


def test_haploid_chrx(fixtures_dir: Path):
    found = False
    with open_variant_file(fixtures_dir / "sample.FB.vcf") as vf:
        for rec in vf:
            if rec.contig == "chrX":
                gt, phased, ploidy = genotype_string(rec, "SAMPLE")
                assert ploidy == 1
                assert gt == "1"
                found = True
                break
    assert found, "chrX record not found"


def test_missing_gt(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "miss.vcf.gz"
    result = harmonize_vcfs([fixtures_dir / "missing_gt.DV.vcf"], out, reference=mini_ref)
    assert result["variant_count"] == 1
