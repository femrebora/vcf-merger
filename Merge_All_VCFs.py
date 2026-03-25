# Merge_All_VCFs.py
#
# Merge a manually specified list of per-caller VCF files from nf-core/sarek
# (FreeBayes, GATK HaplotypeCaller, DeepVariant, and optionally Mutect2 / Strelka)
# into a single output VCF for one patient.
#
# Priority order (highest to lowest): MU > FB > HC > ST > DV
# The file listed first in vcf_files is treated as the highest priority.
#
# Usage:
#   1. Edit the vcf_files list and output_file path below.
#   2. Run:  python Merge_All_VCFs.py
#
# See vcf_utils.py for full merge strategy documentation.

from vcf_utils import merge_vcfs

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Edit the paths below for each patient run.
    # List files in descending priority order (first = highest priority).
    # The internal sort_vcf_files() function will re-sort them by the
    # recognised caller suffixes, so order here only matters for callers
    # that share the same priority level.
    # ------------------------------------------------------------------
    vcf_files = [
        "/path/to/patient/PATIENT_ID.FB.vcf",   # FreeBayes
        "/path/to/patient/PATIENT_ID.HC.vcf",   # GATK HaplotypeCaller
        "/path/to/patient/PATIENT_ID.DV.vcf",   # DeepVariant
    ]

    output_file = "/path/to/output/PATIENT_ID.Merged.vcf"
    # ------------------------------------------------------------------

    merge_vcfs(vcf_files, output_file)
    print(f"Merged VCF saved: {output_file}")
