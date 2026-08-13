"""CLI smoke tests."""

from __future__ import annotations

import json
from pathlib import Path

from vcf_merger.cli import main


def test_cli_inspect(fixtures_dir: Path, capsys):
    rc = main(["inspect", str(fixtures_dir / "sample.HC.vcf")])
    assert rc == 0
    data = json.loads(capsys.readouterr().out)
    assert data["caller"] == "haplotypecaller"


def test_cli_merge(fixtures_dir: Path, mini_ref: Path, tmp_path: Path):
    out = tmp_path / "out.vcf.gz"
    rc = main(
        [
            "merge",
            "--mode",
            "germline",
            "--reference",
            str(mini_ref),
            "-i",
            str(fixtures_dir / "sample.FB.vcf"),
            "-i",
            str(fixtures_dir / "sample.HC.vcf"),
            "-o",
            str(out),
            "--strategy",
            "union",
        ]
    )
    assert rc == 0
    assert out.exists()
