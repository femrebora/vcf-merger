"""Unit tests: canonical identity, gVCF detection, caller detection, validation."""

from __future__ import annotations

from pathlib import Path

import pytest

from vcf_merger.callers import detect_caller, filename_caller_guess, short_label
from vcf_merger.canonical import is_supported_small_variant, make_canonical
from vcf_merger.exceptions import GVcfRejectedError, SampleMismatchError, ValidationError
from vcf_merger.inspect import inspect_vcf
from vcf_merger.models import DetectionSource
from vcf_merger.validation import collect_metadata, reject_gvcf_if_needed, validate_input_set
from vcf_merger.models import AnalysisMode


def test_canonical_identity_case_normalize():
    a = make_canonical("GRCh38", "chr1", 100, "t", "a")
    b = make_canonical("GRCh38", "chr1", 100, "T", "A")
    assert a.key == b.key


def test_symbolic_rejected():
    ok, reason = is_supported_small_variant("T", "<DEL>")
    assert not ok
    assert "symbolic" in (reason or "")


def test_gvcf_detected(fixtures_dir: Path):
    meta = collect_metadata(fixtures_dir / "sample.HC.g.vcf")
    assert meta.is_gvcf
    with pytest.raises(GVcfRejectedError):
        reject_gvcf_if_needed(meta)


def test_caller_detection_header_and_norm_filename(fixtures_dir: Path):
    fb = detect_caller(fixtures_dir / "sample.FB.vcf")
    assert fb.name == "freebayes"
    assert fb.source == DetectionSource.HEADER

    norm = detect_caller(fixtures_dir / "sample.FB.norm.vcf")
    assert norm.name == "freebayes"

    unknown = detect_caller(fixtures_dir / "sample.unknown.vcf")
    assert unknown.name == "unknown"

    explicit = detect_caller(fixtures_dir / "sample.unknown.vcf", explicit="DV")
    assert explicit.name == "deepvariant"
    assert explicit.source == DetectionSource.CLI


def test_filename_guess_not_broken_by_norm_suffix(fixtures_dir: Path):
    assert filename_caller_guess(fixtures_dir / "sample.FB.norm.vcf") == "freebayes"


def test_assembly_mismatch(fixtures_dir: Path, mini_ref: Path):
    metas = [
        inspect_vcf(fixtures_dir / "sample.FB.vcf"),
        inspect_vcf(fixtures_dir / "hg19.FB.vcf"),
    ]
    with pytest.raises(ValidationError, match="assembl"):
        validate_input_set(metas, mode=AnalysisMode.GERMLINE, reference=str(mini_ref))


def test_sample_mismatch(fixtures_dir: Path, mini_ref: Path):
    metas = [
        inspect_vcf(fixtures_dir / "family.HC.vcf"),
        inspect_vcf(fixtures_dir / "family_mismatch.FB.vcf"),
    ]
    with pytest.raises(SampleMismatchError):
        validate_input_set(metas, mode=AnalysisMode.GERMLINE, reference=str(mini_ref))


def test_contig_length_mismatch(fixtures_dir: Path, mini_ref: Path):
    metas = [
        inspect_vcf(fixtures_dir / "sample.FB.vcf"),
        inspect_vcf(fixtures_dir / "bad_contig.HC.vcf"),
    ]
    with pytest.raises(ValidationError, match="contig"):
        validate_input_set(metas, mode=AnalysisMode.GERMLINE, reference=str(mini_ref))


def test_inspect_gz(fixtures_dir: Path):
    meta = inspect_vcf(fixtures_dir / "sample.FB.vcf.gz")
    assert meta.caller == "freebayes"
    assert "SAMPLE" in meta.samples


def test_short_label():
    assert short_label("freebayes") == "FB"
    assert short_label("haplotypecaller") == "HC"
