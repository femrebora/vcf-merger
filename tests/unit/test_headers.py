"""Header reconciliation tests."""

from __future__ import annotations

from pathlib import Path

import pysam

from vcf_merger.headers import reconcile_headers
from vcf_merger.io_vcf import open_variant_file


def test_header_conflict_namespaced(fixtures_dir: Path):
    with open_variant_file(fixtures_dir / "hdr_conflict_a.vcf") as a, open_variant_file(
        fixtures_dir / "hdr_conflict_b.vcf"
    ) as b:
        out = reconcile_headers([a.header, b.header], samples=["SAMPLE"])
    info_ids = set(out.info.keys())
    assert "XX" in info_ids
    assert any(i.startswith("XX_VM") for i in info_ids)
    assert "VM_CALLERS" in info_ids
    assert "VM_GT_CONFLICT" in info_ids
