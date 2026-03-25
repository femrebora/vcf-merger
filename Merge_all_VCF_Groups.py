# Merge_all_VCF_Groups.py
#
# Batch-merge per-caller VCF files for all patients in a directory.
# Files are grouped by patient ID — the part of the filename before the first dot.
#
# Example layout:
#   PATIENT_01.FB.vcf  ──┐
#   PATIENT_01.HC.vcf  ──┤──► PATIENT_01.Merged.vcf
#   PATIENT_01.DV.vcf  ──┘
#
# Groups with 2 or fewer files are skipped (a patient must have output from
# at least 3 callers to be processed).
#
# Priority order (highest to lowest): MU > FB > HC > ST > DV
#
# Usage:
#   python Merge_all_VCF_Groups.py <input_directory> <output_directory>
#
# Example:
#   python Merge_all_VCF_Groups.py /data/sarek_output /data/merged_vcfs
#
# See vcf_utils.py for full merge strategy documentation.

import os
import argparse
from glob import glob

from vcf_utils import merge_vcfs


def find_vcf_groups(directory: str) -> dict[str, list[str]]:
    """Group VCF files by patient ID (filename prefix before the first dot)."""
    vcf_files = glob(os.path.join(directory, '*.vcf'))
    vcf_groups: dict[str, list[str]] = {}

    for vcf in vcf_files:
        base_name = os.path.basename(vcf)
        group_key = base_name.split('.')[0]
        if group_key not in vcf_groups:
            vcf_groups[group_key] = []
        vcf_groups[group_key].append(vcf)

    return vcf_groups


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Batch-merge VCF files by patient group from an nf-core/sarek run.'
    )
    parser.add_argument('directory', type=str, help='Directory containing input VCF files')
    parser.add_argument('output_dir', type=str, help='Directory to save merged VCF files')
    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    vcf_groups = find_vcf_groups(args.directory)

    for group, files in vcf_groups.items():
        if len(files) > 2:
            output_file = os.path.join(args.output_dir, f"{group}.Merged.vcf")
            merge_vcfs(files, output_file)
            print(f"Merged VCF saved: {output_file}")
        else:
            print(
                f"Skipping '{group}': only {len(files)} VCF file(s) found "
                f"(expected more than 2)."
            )
