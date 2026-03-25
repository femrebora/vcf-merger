# vcf-merger

Python scripts for merging per-caller VCF files produced by the [nf-core/sarek](https://nf-co.re/sarek) germline variant calling pipeline.

## Background

nf-core/sarek runs multiple variant callers in parallel on the same patient sample:

| Caller | File suffix | Description |
|--------|-------------|-------------|
| [FreeBayes](https://github.com/freebayes/freebayes) | `.FB.vcf` | Bayesian haplotype-based caller |
| [GATK HaplotypeCaller](https://gatk.broadinstitute.org/hc/en-us/articles/360037225632) | `.HC.vcf` | Graph-based local assembly caller |
| [DeepVariant](https://github.com/google/deepvariant) | `.DV.vcf` | Deep-learning-based caller |
| [Mutect2](https://gatk.broadinstitute.org/hc/en-us/articles/360037593851) *(optional)* | `.MU.vcf` | Somatic/germline caller |
| [Strelka](https://github.com/Illumina/strelka) *(optional)* | `.ST.vcf` | Fast germline/somatic caller |

Each caller has different strengths. These scripts merge their outputs into a single VCF per patient, using a **priority-based strategy**: when the same variant is called by multiple tools, metadata columns (QUAL, FILTER, INFO, FORMAT, genotype) from the higher-priority caller are kept. Lower-priority callers only fill in values that are missing (`NaN`) in the higher-priority result.

**Default priority order:** MU > FB > HC > ST > DV

Variants are matched on `CHROM + POS + REF + ALT`, so multi-allelic sites (same genomic position, different alleles) are always kept as separate rows.

---

## Scripts

### `Merge_All_VCFs.py` — single patient, manual paths

Use this when you want to merge VCF files for **one patient** by explicitly listing the file paths.

**Edit the file:**
```python
vcf_files = [
    "/path/to/PATIENT_ID.FB.vcf",   # highest priority
    "/path/to/PATIENT_ID.HC.vcf",
    "/path/to/PATIENT_ID.DV.vcf",   # lowest priority
]
output_file = "/path/to/PATIENT_ID.Merged.vcf"
```

**Run:**
```bash
python Merge_All_VCFs.py
```

---

### `Merge_all_VCF_Groups.py` — batch, whole directory

Use this when you have a directory containing VCF files for **multiple patients** and want to process them all at once.

Files are grouped by patient ID — the part of the filename before the first dot:

```
PATIENT_01.FB.vcf  ──┐
PATIENT_01.HC.vcf  ──┤──► PATIENT_01.Merged.vcf
PATIENT_01.DV.vcf  ──┘

PATIENT_02.FB.vcf  ──┐
PATIENT_02.HC.vcf  ──┤──► PATIENT_02.Merged.vcf
PATIENT_02.DV.vcf  ──┘
```

Groups with 2 or fewer files are skipped with a warning (a patient must have output from at least 3 callers).

**Run:**
```bash
python Merge_all_VCF_Groups.py /data/sarek_output /data/merged_vcfs
```

| Argument | Description |
|----------|-------------|
| `directory` | Directory containing the input `.vcf` files |
| `output_dir` | Directory where merged VCF files will be written (created if it does not exist) |

---

## How the merge works

```
For each patient group:
  1. Sort files by caller priority (MU > FB > HC > ST > DV)
  2. Start with the highest-priority caller as the base DataFrame
  3. For each subsequent caller:
       - Outer join on [CHROM, POS, REF, ALT]
         → variants unique to either caller are included (no variants are lost)
       - For shared columns (QUAL, FILTER, INFO, FORMAT, genotype …):
           keep the higher-priority value; fill NaN from the lower-priority value
  4. Drop exact duplicates on [CHROM, POS, ID, REF, ALT]
  5. Write result as tab-separated file
```

---

## Installation

```bash
pip install -r requirements.txt
```

Requires Python 3.8+.

---

## Known limitations

- **VCF meta-information lines (`##`) are not written to the output.** The output file contains only the column header line (`#CHROM POS ...`) and data rows. If your downstream tools require a fully compliant VCF header, add the `##` lines manually or use a tool like `bcftools reheader`.
- Caller-specific FORMAT fields (e.g. DeepVariant's `VAF` tag vs HaplotypeCaller's `AD` tag) may differ. The merge preserves the FORMAT string and sample column from the highest-priority caller that called the variant. Fields unique to lower-priority callers are not carried over.
- This tool does not perform genotype-level reconciliation. For clinical-grade ensemble calling, consider [GATK CombineGVCFs](https://gatk.broadinstitute.org/hc/en-us/articles/360037053272) or [GLnexus](https://github.com/dnanexus-rnd/GLnexus).

---

## License

MIT
