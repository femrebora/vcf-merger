"""Unit tests for genotype reconciliation and ensemble strategy."""

from __future__ import annotations

from vcf_merger.harmonizer import passes_strategy, reconcile_genotypes
from vcf_merger.models import CallerEvidence, EnsembleStrategy, SampleEvidence


def _ev(caller: str, gt: str, passed: bool = True, sample: str = "SAMPLE") -> CallerEvidence:
    return CallerEvidence(
        caller=caller,
        caller_version=None,
        source_file="x.vcf",
        filter_status="PASS" if passed else "LowQual",
        passed=passed,
        qual=10.0,
        samples=[SampleEvidence(sample=sample, genotype=gt, ploidy=2)],
    )


def test_reconcile_unanimous():
    gt, conflict = reconcile_genotypes(
        [_ev("freebayes", "0/1"), _ev("haplotypecaller", "0|1")],
        "SAMPLE",
    )
    assert conflict is False
    assert gt in ("0/1", "0|1")


def test_reconcile_conflict():
    gt, conflict = reconcile_genotypes(
        [_ev("freebayes", "0/1"), _ev("haplotypecaller", "1/1")],
        "SAMPLE",
    )
    assert conflict is True
    assert gt is not None


def test_strategy_union_includes_filtered():
    ev = [_ev("freebayes", "0/1", passed=False)]
    assert passes_strategy(ev, EnsembleStrategy.UNION)
    assert not passes_strategy(ev, EnsembleStrategy.PASS_UNION)


def test_strategy_consensus():
    ev = [_ev("freebayes", "0/1"), _ev("deepvariant", "0/1")]
    assert passes_strategy(ev, EnsembleStrategy.CONSENSUS, consensus_n=2)
    assert not passes_strategy(ev[:1], EnsembleStrategy.CONSENSUS, consensus_n=2)
