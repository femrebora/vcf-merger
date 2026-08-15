# Merge_all_VCF_Groups.py
#
# Legacy batch wrapper around vcf_merger.harmonize_vcfs.
# Groups files by the filename prefix before the first dot.
# The old "require >2 callers" gate has been removed.
#
# Usage:
#   python Merge_all_VCF_Groups.py <input_directory> <output_directory> <reference.fasta>

from __future__ import annotations

import argparse
import os
import warnings
from glob import glob

from vcf_merger.harmonizer import harmonize_vcfs

warnings.warn(
    "Merge_all_VCF_Groups.py is deprecated; use the vcf-merger CLI. See docs/migration.md.",
    DeprecationWarning,
    stacklevel=1,
)


def find_vcf_groups(directory: str) -> dict[str, list[str]]:
    patterns = ["*.vcf", "*.vcf.gz"]
    vcf_files: list[str] = []
    for pat in patterns:
        vcf_files.extend(glob(os.path.join(directory, pat)))
    # Skip secondary gVCF naming in default discovery by leaving them to validation
    vcf_groups: dict[str, list[str]] = {}
    for vcf in vcf_files:
        base_name = os.path.basename(vcf)
        if ".g.vcf" in base_name.lower():
            continue
        group_key = base_name.split(".")[0]
        vcf_groups.setdefault(group_key, []).append(vcf)
    return vcf_groups


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Batch-harmonize per-caller VCF files grouped by sample ID prefix."
    )
    parser.add_argument("directory", type=str, help="Directory containing input VCF files")
    parser.add_argument("output_dir", type=str, help="Directory for harmonized outputs")
    parser.add_argument("reference", type=str, help="Reference FASTA for normalization")
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    vcf_groups = find_vcf_groups(args.directory)

    for group, files in vcf_groups.items():
        if len(files) < 1:
            continue
        output_file = os.path.join(args.output_dir, f"{group}.harmonized.vcf.gz")
        result = harmonize_vcfs(
            files,
            output_file,
            reference=args.reference,
            mode="germline",
            strategy="union",
            command_line=["Merge_all_VCF_Groups.py", args.directory, args.output_dir, args.reference],
        )
        print(f"Merged VCF saved: {result['paths']['vcf']}")
