"""Lightweight performance / scaling check."""

from __future__ import annotations

import resource
import time
from pathlib import Path

import pytest

from vcf_merger.harmonizer import harmonize_vcfs

pytestmark = pytest.mark.performance


def _write_large_vcf(path: Path, n: int, caller_tag: str) -> None:
    with path.open("w") as fh:
        fh.write("##fileformat=VCFv4.2\n")
        fh.write("##assembly=GRCh38\n")
        fh.write(f"##source={caller_tag}\n")
        fh.write("##contig=<ID=chr1,length=2000000>\n")
        fh.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        fh.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="DP">\n')
        fh.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        # Use positions that exist on mini_ref pattern - but we'll use --no-normalize
        # with a fake ref-less merge for throughput of the merge loop.
        for i in range(n):
            pos = 100 + i
            fh.write(f"chr1\t{pos}\t.\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30\n")


@pytest.mark.performance
def test_merge_throughput_three_callers(tmp_path: Path):
    n = 5000
    paths = []
    for tag in ("freebayes", "HaplotypeCaller", "DeepVariant"):
        p = tmp_path / f"s.{tag}.vcf"
        _write_large_vcf(p, n, tag)
        paths.append(p)
    out = tmp_path / "out.vcf.gz"
    rss_before = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    t0 = time.perf_counter()
    result = harmonize_vcfs(
        paths,
        out,
        normalize=False,
        strategy="union",
        mode="germline",
    )
    elapsed = time.perf_counter() - t0
    rss_after = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    assert result["variant_count"] == n
    rps = n / elapsed if elapsed else float("inf")
    # Soft assertion: should process at a reasonable rate in CI
    assert rps > 100
    # Document memory delta (Linux ru_maxrss is KB)
    mem_delta_kb = max(0, rss_after - rss_before)
    (tmp_path / "perf.txt").write_text(
        f"records={n} callers=3 seconds={elapsed:.3f} rps={rps:.1f} maxrss_delta_kb={mem_delta_kb}\n"
    )
