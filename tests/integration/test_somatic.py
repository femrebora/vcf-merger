"""Somatic-mode unit and integration tests."""

from __future__ import annotations

import json
from pathlib import Path

import pysam
import pytest

from vcf_merger.af_utils import strelka_snv_af
from vcf_merger.cli import main
from vcf_merger.exceptions import ValidationError
from vcf_merger.harmonizer import harmonize_vcfs
from vcf_merger.somatic import resolve_somatic_context
from vcf_merger.inspect import inspect_vcf


def test_strelka_snv_af_helper():
    fmt = {"AU": [10, 12], "CU": [0, 0], "GU": [0, 0], "TU": [40, 42]}
    assert strelka_snv_af(fmt, "A") == pytest.approx(0.2)


def test_resolve_tumor_normal(fixtures_dir: Path):
    metas = [
        inspect_vcf(fixtures_dir / "somatic.MU.vcf"),
        inspect_vcf(fixtures_dir / "somatic.ST.vcf"),
    ]
    ctx = resolve_somatic_context(metas, tumor_sample="TUMOR", normal_sample="NORMAL")
    assert ctx.tumor_sample == "TUMOR"
    assert ctx.normal_sample == "NORMAL"
    assert ctx.tumor_only is False


def test_resolve_tumor_only(fixtures_dir: Path):
    metas = [inspect_vcf(fixtures_dir / "somatic_tumor_only.MU.vcf")]
    ctx = resolve_somatic_context(metas, tumor_sample="TUMOR", tumor_only=True)
    assert ctx.tumor_only is True
    assert ctx.normal_sample is None


def test_tumor_only_rejects_normal(fixtures_dir: Path):
    metas = [inspect_vcf(fixtures_dir / "somatic_tumor_only.MU.vcf")]
    with pytest.raises(ValidationError):
        resolve_somatic_context(
            metas, tumor_sample="TUMOR", normal_sample="NORMAL", tumor_only=True
        )


@pytest.mark.bcftools
def test_somatic_merge_afs(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "som.vcf.gz"
    result = harmonize_vcfs(
        [fixtures_dir / "somatic.MU.vcf", fixtures_dir / "somatic.ST.vcf"],
        out,
        reference=mini_ref,
        mode="somatic",
        tumor_sample="TUMOR",
        normal_sample="NORMAL",
        strategy="union",
    )
    assert result["somatic"]["tumor_sample"] == "TUMOR"
    assert result["somatic"]["normal_sample"] == "NORMAL"

    with pysam.VariantFile(out) as vf:
        recs = { (r.pos, r.alts[0]): r for r in vf }
    hit = recs[(100, "A")]
    assert "VM_SOMATIC_PAIR" in hit.info
    assert hit.info["VM_SOMATIC_PAIR"] == "TUMOR:NORMAL"
    assert "VM_TUMOR_AF" in hit.info
    assert float(hit.info["VM_TUMOR_AF"]) > 0
    assert "VM_NORMAL_AF" in hit.info

    with open(result["paths"]["evidence"], encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh]
    row = next(r for r in rows if r["canonical"]["pos"] == 100)
    assert row["somatic"]["tumor_sample"] == "TUMOR"
    st = next(e for e in row["evidence"] if e["caller"] == "strelka")
    assert st["tumor_allele_fraction"] == pytest.approx(0.2)
    mu = next(e for e in row["evidence"] if e["caller"] == "mutect2")
    assert mu["normal_allele_fraction"] == pytest.approx(0.0)


@pytest.mark.bcftools
def test_somatic_tumor_only_merge(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "to.vcf.gz"
    result = harmonize_vcfs(
        [fixtures_dir / "somatic_tumor_only.MU.vcf"],
        out,
        reference=mini_ref,
        mode="somatic",
        tumor_sample="TUMOR",
        tumor_only=True,
        strategy="pass-union",
    )
    assert result["somatic"]["tumor_only"] is True
    assert result["variant_count"] == 2
    with pysam.VariantFile(out) as vf:
        for r in vf:
            assert r.info.get("VM_SOMATIC_PAIR") == "TUMOR"
            assert "VM_TUMOR_AF" in r.info


@pytest.mark.bcftools
def test_cli_somatic(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "cli_som.vcf.gz"
    ev = tmp_path / "custom.evidence.jsonl"
    prov = tmp_path / "custom.provenance.json"
    rc = main(
        [
            "merge",
            "--mode",
            "somatic",
            "--reference",
            str(mini_ref),
            "-i",
            str(fixtures_dir / "somatic.MU.vcf"),
            "-i",
            str(fixtures_dir / "somatic.ST.vcf"),
            "--tumor-sample",
            "TUMOR",
            "--normal-sample",
            "NORMAL",
            "-o",
            str(out),
            "--evidence-output",
            str(ev),
            "--provenance-output",
            str(prov),
        ]
    )
    assert rc == 0
    assert out.exists() and ev.exists() and prov.exists()
    data = json.loads(prov.read_text())
    assert data["analysis_mode"] == "somatic"
    assert data["normalization"]["somatic"]["tumor_sample"] == "TUMOR"
