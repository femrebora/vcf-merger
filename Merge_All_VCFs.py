# Merge_All_VCFs.py
#
# Legacy wrapper around vcf_merger.harmonize_vcfs.
# Prefer: vcf-merger merge --mode germline --reference REF.fa -i ... -o ...
#
# Edit the paths below, then run:  python Merge_All_VCFs.py

from __future__ import annotations

import warnings

from vcf_merger.harmonizer import harmonize_vcfs

warnings.warn(
    "Merge_All_VCFs.py is deprecated; use the vcf-merger CLI or vcf_merger.harmonize_vcfs. "
    "See docs/migration.md.",
    DeprecationWarning,
    stacklevel=1,
)

if __name__ == "__main__":
    # ------------------------------------------------------------------
    # Edit paths for each run. Provide a reference FASTA for normalization.
    # ------------------------------------------------------------------
    vcf_files = [
        "/path/to/patient/PATIENT_ID.FB.vcf",
        "/path/to/patient/PATIENT_ID.HC.vcf",
        "/path/to/patient/PATIENT_ID.DV.vcf",
    ]
    output_file = "/path/to/output/PATIENT_ID.harmonized.vcf.gz"
    reference = "/path/to/reference.fasta"
    # ------------------------------------------------------------------

    result = harmonize_vcfs(
        vcf_files,
        output_file,
        reference=reference,
        mode="germline",
        strategy="union",
        command_line=["Merge_All_VCFs.py"],
    )
    print(f"Merged VCF saved: {result['paths']['vcf']}")
    print(f"Evidence: {result['paths']['evidence']}")
    print(f"Provenance: {result['paths']['provenance']}")
