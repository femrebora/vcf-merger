# Merge_all_VCF_Groups
#
# Batch-processes a directory of VCF files produced by different variant callers
# (FreeBayes .FB.vcf, GATK HaplotypeCaller .HC.vcf, DeepVariant .DV.vcf, and
# optionally Mutect2 .MU.vcf and Strelka .ST.vcf) from the nf-core/sarek pipeline.
#
# Files are grouped by patient ID (the part of the filename before the first dot).
# Example: PATIENT_01.FB.vcf, PATIENT_01.HC.vcf, PATIENT_01.DV.vcf
#          → group key: PATIENT_01
#
# Priority order (highest to lowest):
#   MU (Mutect2) > FB (FreeBayes) > HC (HaplotypeCaller) > ST (Strelka) > DV (DeepVariant)
#
# A group must contain more than 2 files to be merged; groups with 1 or 2 files are
# skipped with a warning. This guard ensures that only patients with calls from all
# expected callers are processed.
#
# Usage:
#   python Merge_all_VCF_Groups.py <input_directory> <output_directory>
#
# Example:
#   python Merge_all_VCF_Groups.py /data/sarek_output /data/merged_vcfs
#
# Note: VCF meta-information lines (##) are stripped. The output file contains the
#       column header (#CHROM ...) and data rows only.

import os
import io
import argparse
from glob import glob

import pandas as pd


def read_vcf(file):
    """Read a VCF file into a pandas DataFrame, skipping ## meta-information lines."""
    with open(file, 'r', encoding='iso-8859-1') as f:
        lines = [line for line in f if not line.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t',
        dtype=str,        # keep all values as strings to avoid type coercion
        low_memory=False
    )
    return df


def sort_vcf_files(vcf_files):
    """Sort VCF files by caller priority (MU > FB > HC > ST > DV)."""
    priority_order = ['.MU.vcf', '.FB.vcf', '.HC.vcf', '.ST.vcf', '.DV.vcf']
    vcf_files.sort(
        key=lambda x: next(
            (i for i, ext in enumerate(priority_order) if ext in x),
            len(priority_order)   # unknown callers go last
        )
    )
    return vcf_files


def merge_vcfs(vcf_files, output_file):
    """Merge multiple VCF files with caller priority.

    Variants are matched on CHROM + POS + REF + ALT so that multi-allelic sites
    (same position, different alleles) are kept as separate rows instead of being
    incorrectly collapsed.
    """
    vcf_files = sort_vcf_files(vcf_files)
    dfs = [read_vcf(vcf) for vcf in vcf_files]

    merged_df = dfs[0]
    for df in dfs[1:]:
        merged_df = pd.merge(
            merged_df, df,
            on=['#CHROM', 'POS', 'REF', 'ALT'],
            how='outer',
            suffixes=('', '_y')
        )
        # For every column that duplicated during the merge, keep the value from
        # the higher-priority (left) dataframe and fall back to the lower-priority
        # (right) value only when the left is missing.
        for col in list(merged_df.columns):
            if col.endswith('_y'):
                original_col = col[:-2]          # strip the '_y' suffix
                merged_df[original_col] = merged_df[original_col].combine_first(
                    merged_df[col]
                )
                merged_df.drop(columns=[col], inplace=True)

    # Safety deduplication: drop any remaining identical variant rows
    merged_df = merged_df.drop_duplicates(subset=['#CHROM', 'POS', 'ID', 'REF', 'ALT'])

    merged_df.to_csv(output_file, sep='\t', index=False)


def find_vcf_groups(directory):
    """Group VCF files by patient ID (filename prefix before the first dot)."""
    vcf_files = glob(os.path.join(directory, '*.vcf'))
    vcf_groups = {}

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
