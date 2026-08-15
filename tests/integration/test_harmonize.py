"""Integration tests for harmonize / normalize (require bcftools)."""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pysam
import pytest

from vcf_merger.exceptions import GVcfRejectedError, ValidationError
from vcf_merger.harmonizer import harmonize_vcfs
from vcf_merger.normalization import normalize_vcf

bcftools = shutil.which("bcftools")
pytestmark = pytest.mark.skipif(bcftools is None, reason="bcftools not available")


@pytest.mark.bcftools
def test_normalize_idempotent(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out1 = tmp_path / "n1.vcf.gz"
    out2 = tmp_path / "n2.vcf.gz"
    normalize_vcf(fixtures_dir / "sample.FB.vcf", out1, reference=mini_ref)
    normalize_vcf(out1, out2, reference=mini_ref)
    with pysam.VariantFile(out1) as a, pysam.VariantFile(out2) as b:
        ra = [(r.contig, r.pos, r.ref, r.alts) for r in a]
        rb = [(r.contig, r.pos, r.ref, r.alts) for r in b]
    assert ra == rb


@pytest.mark.bcftools
def test_merge_union_preserves_evidence_order_independent(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out_a = tmp_path / "a.harmonized.vcf.gz"
    out_b = tmp_path / "b.harmonized.vcf.gz"
    inputs_a = [
        fixtures_dir / "sample.FB.vcf",
        fixtures_dir / "sample.HC.vcf",
        fixtures_dir / "sample.DV.vcf",
    ]
    inputs_b = list(reversed(inputs_a))
    ra = harmonize_vcfs(inputs_a, out_a, reference=mini_ref, mode="germline", strategy="union")
    rb = harmonize_vcfs(inputs_b, out_b, reference=mini_ref, mode="germline", strategy="union")
    assert ra["variant_count"] == rb["variant_count"]

    def keys(path: Path):
        with pysam.VariantFile(path) as vf:
            return sorted(
                (
                    r.contig,
                    r.pos,
                    r.ref,
                    r.alts[0],
                    r.info.get("VM_CALLER_COUNT"),
                    tuple(sorted(str(r.info.get("VM_CALLERS", "")).split(","))),
                )
                for r in vf
            )

    assert keys(out_a) == keys(out_b)

    # Evidence sidecar has all three callers for chr1:100 T>A
    with open(ra["paths"]["evidence"], encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh]
    hit = [
        r
        for r in rows
        if r["canonical"]["pos"] == 100 and r["canonical"]["ref"] == "T" and r["canonical"]["alt"] == "A"
    ]
    assert hit
    assert set(hit[0]["callers_all"]) >= {"freebayes", "haplotypecaller", "deepvariant"}


@pytest.mark.bcftools
def test_gt_conflict_flag(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "conflict.vcf.gz"
    harmonize_vcfs(
        [fixtures_dir / "conflict.FB.vcf", fixtures_dir / "conflict.HC.vcf"],
        out,
        reference=mini_ref,
        strategy="union",
    )
    with pysam.VariantFile(out) as vf:
        recs = list(vf)
    assert len(recs) == 1
    assert "VM_GT_CONFLICT" in recs[0].info


@pytest.mark.bcftools
def test_pass_union_filters(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    # Only FB LowQual indel at 300 should drop under pass-union if no other caller has it
    out = tmp_path / "pass.vcf.gz"
    result = harmonize_vcfs(
        [fixtures_dir / "sample.FB.vcf", fixtures_dir / "sample.HC.vcf"],
        out,
        reference=mini_ref,
        strategy="pass-union",
    )
    with pysam.VariantFile(out) as vf:
        positions = {(r.contig, r.pos, r.alts[0]) for r in vf}
    # chr1:300 deletion is LowQual in FB only — excluded
    assert not any(p[1] == 300 for p in positions)
    assert result["variant_count"] >= 1


@pytest.mark.bcftools
def test_gvcf_rejected_on_merge(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    with pytest.raises(GVcfRejectedError):
        harmonize_vcfs(
            [fixtures_dir / "sample.HC.g.vcf", fixtures_dir / "sample.FB.vcf"],
            tmp_path / "out.vcf.gz",
            reference=mini_ref,
        )


@pytest.mark.bcftools
def test_sv_not_in_vcf_but_in_sidecar(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "sv.vcf.gz"
    result = harmonize_vcfs(
        [fixtures_dir / "sv.FB.vcf"],
        out,
        reference=mini_ref,
        strategy="union",
    )
    with pysam.VariantFile(out) as vf:
        alts = [r.alts[0] for r in vf]
    assert "<DEL>" not in alts
    assert "G" in alts
    with open(result["paths"]["evidence"], encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh]
    assert any(r.get("unsupported") for r in rows)


@pytest.mark.bcftools
def test_empty_vcf(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "empty.out.vcf.gz"
    result = harmonize_vcfs([fixtures_dir / "empty.vcf"], out, reference=mini_ref)
    assert result["variant_count"] == 0
    assert Path(result["paths"]["provenance"]).exists()


@pytest.mark.bcftools
def test_multiallelic_split(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "multi.vcf.gz"
    harmonize_vcfs([fixtures_dir / "multiallelic.FB.vcf"], out, reference=mini_ref)
    with pysam.VariantFile(out) as vf:
        alts = sorted(r.alts[0] for r in vf)
    assert alts == ["A", "G"]


@pytest.mark.bcftools
def test_roundtrip_pysam(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "rt.vcf.gz"
    harmonize_vcfs(
        [fixtures_dir / "sample.FB.vcf", fixtures_dir / "sample.HC.vcf"],
        out,
        reference=mini_ref,
    )
    with pysam.VariantFile(out) as vf:
        assert vf.header is not None
        list(vf)  # round-trip read
    assert (tmp_path / "rt.vcf.gz.tbi").exists() or (tmp_path / "rt.vcf.gz.csi").exists()


@pytest.mark.bcftools
def test_somatic_mode(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
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
    assert result["mode"] == "somatic"
    with open(result["paths"]["evidence"], encoding="utf-8") as fh:
        rows = [json.loads(line) for line in fh]
    hit = [r for r in rows if r["canonical"]["pos"] == 100]
    assert hit
    callers = set(hit[0]["callers_all"])
    assert "mutect2" in callers and "strelka" in callers
    # Mutect2 metrics preserved in sidecar
    mu = next(e for e in hit[0]["evidence"] if e["caller"] == "mutect2")
    assert "TLOD" in mu["caller_specific_metrics"]


@pytest.mark.bcftools
def test_bcftools_view_accepts_output(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    import subprocess

    out = tmp_path / "valid.vcf.gz"
    harmonize_vcfs(
        [fixtures_dir / "sample.FB.vcf", fixtures_dir / "sample.DV.vcf"],
        out,
        reference=mini_ref,
    )
    proc = subprocess.run(["bcftools", "view", "-H", str(out)], capture_output=True, text=True)
    assert proc.returncode == 0
