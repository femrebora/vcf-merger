# Merge All VCF Files
#
# Merges a manually specified list of VCF files produced by different variant callers
# (e.g. FreeBayes, GATK HaplotypeCaller, DeepVariant) from the nf-core/sarek pipeline.
#
# Priority: the first file in vcf_files list has the highest priority. When the same
# variant (same CHROM, POS, REF, ALT) is called by multiple tools, columns from the
# higher-priority file are kept; columns from lower-priority files only fill in where
# the higher-priority file has missing (NaN) values.
#
# Usage:
#   1. Edit the vcf_files list below with absolute paths to your patient's VCF files.
#   2. Edit output_file with the desired output path.
#   3. Run:  python Merge_All_VCFs.py
#
# Note: VCF meta-information lines (##) are stripped. The output file contains the
#       column header (#CHROM ...) and data rows only.

import pandas as pd
import io


def read_vcf(file):
    """Read a VCF file into a pandas DataFrame, skipping ## meta-information lines."""
    with open(file, 'r') as f:
        lines = [line for line in f if not line.startswith('##')]
    df = pd.read_csv(
        io.StringIO(''.join(lines)),
        sep='\t'
    )
    return df


def merge_vcfs(vcf_files, output_file):
    """Merge multiple VCF files with caller priority.

    Variants are matched on CHROM + POS + REF + ALT so that multi-allelic sites
    (same position, different alleles) are kept as separate rows instead of being
    incorrectly collapsed.
    """
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

    merged_df.to_csv(output_file, sep='\t', index=False)


if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Edit the paths below for each patient run.
    # List files in priority order: first file = highest priority.
    # ------------------------------------------------------------------
    vcf_files = [
        "/path/to/patient/PATIENT_ID.FB.vcf",   # FreeBayes (highest priority)
        "/path/to/patient/PATIENT_ID.HC.vcf",   # GATK HaplotypeCaller
        "/path/to/patient/PATIENT_ID.DV.vcf",   # DeepVariant (lowest priority)
    ]

    output_file = "/path/to/output/PATIENT_ID.Merged.vcf"
    # ------------------------------------------------------------------

    merge_vcfs(vcf_files, output_file)
    print(f"Merged VCF file saved as {output_file}")
