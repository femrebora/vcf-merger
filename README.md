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

Each caller has different strengths. These scripts merge their outputs into a single VCF per patient, using a **priority-based strategy**.

**Default priority order:** MU > FB > HC > ST > DV

---

## Merge strategy

For each unique variant identified by `CHROM + POS + REF + ALT`:

1. The **complete row** (QUAL, FILTER, INFO, FORMAT, sample genotype) from the **highest-priority caller** that found the variant is used. This avoids mixing FORMAT tags from different callers, which would produce an internally inconsistent genotype column.
2. A `CALLERS=FB,HC` tag is appended to the INFO field so you always know which callers supported each variant.
3. Variants unique to any single caller are included — **no variants are discarded**.
4. The `##` meta-information header lines from all callers are merged and written at the top of the output, producing a valid VCF file.
5. Output is sorted by chromosome and genomic position.

```
Example output row (variant found by all three callers):
#CHROM  POS     ID  REF  ALT  QUAL  FILTER  INFO                        FORMAT  SAMPLE
chr1    925952  .   G    A    312   PASS    DP=45;AF=0.48;CALLERS=FB,HC,DV  GT:AD  0/1:23,22
                                            ↑ INFO from FB (highest priority)
                                                              ↑ provenance tag
```

---

## Files

| File | Purpose |
|------|---------|
| `vcf_utils.py` | Core merge logic shared by both scripts |
| `Merge_All_VCFs.py` | Merge VCF files for one patient (manually specify paths) |
| `Merge_all_VCF_Groups.py` | Batch merge for all patients in a directory |

---

## Installation

```bash
pip install -r requirements.txt
```

Requires Python 3.10+.

---

## Usage

### Single patient — `Merge_All_VCFs.py`

Edit the file paths at the bottom of `Merge_All_VCFs.py`:

```python
vcf_files = [
    "/path/to/PATIENT_ID.FB.vcf",
    "/path/to/PATIENT_ID.HC.vcf",
    "/path/to/PATIENT_ID.DV.vcf",
]
output_file = "/path/to/PATIENT_ID.Merged.vcf"
```

Then run:

```bash
python Merge_All_VCFs.py
```

---

### Batch — `Merge_all_VCF_Groups.py`

Place all patient VCF files in one directory with the naming convention `PATIENT_ID.CALLER.vcf`:

```
/data/sarek_output/
├── PATIENT_01.FB.vcf
├── PATIENT_01.HC.vcf
├── PATIENT_01.DV.vcf
├── PATIENT_02.FB.vcf
├── PATIENT_02.HC.vcf
└── PATIENT_02.DV.vcf
```

Run:

```bash
python Merge_all_VCF_Groups.py /data/sarek_output /data/merged_vcfs
```

Output:

```
/data/merged_vcfs/
├── PATIENT_01.Merged.vcf
└── PATIENT_02.Merged.vcf
```

Groups with 2 or fewer VCF files are skipped with a warning (a patient must have output from at least 3 callers).

---

## Known limitations

### Summary

| Limitation | Severity | Workaround |
|------------|----------|------------|
| Indel representation differences | High | Run `bcftools norm` on each input VCF before merging |
| CALLERS concordance signal is not equal across callers | Medium | Treat HC-only calls separately (see below) |
| QUAL scores are not comparable between callers | Medium | Use `CALLERS` count and `FILTER=PASS` as quality proxy instead |
| Caller-specific FORMAT fields not carried across callers | Low | FORMAT is always internally consistent per row |
| Not tested on all possible VCF edge cases | Low | Validate output on one patient before batch processing |

---

### 1. Indel representation

Different callers may represent the same indel differently (e.g. different left-alignment or REF/ALT padding). The same true indel can appear as two separate rows if callers disagree on representation. Normalize your input VCFs with `bcftools norm` before running the merge:

```bash
bcftools norm -f reference.fasta -m -any patient.FB.vcf -o patient.FB.norm.vcf
bcftools norm -f reference.fasta -m -any patient.HC.vcf -o patient.HC.norm.vcf
bcftools norm -f reference.fasta -m -any patient.DV.vcf -o patient.DV.norm.vcf
```

---

### 2. CALLERS concordance signal is not equal across callers

The `CALLERS=` tag counts how many independent callers found a variant. However, not all callers are run in the same mode:

- **FreeBayes** and **DeepVariant** are always run independently per sample
- **HaplotypeCaller** may be run in joint genotyping mode, where it sees all samples in the cohort simultaneously and can rescue variants with weak per-sample evidence using cohort-level statistics

This means:

| CALLERS value | Interpretation |
|---------------|---------------|
| `CALLERS=FB,HC,DV` | All three agree — high confidence |
| `CALLERS=FB,DV` | Two independent callers agree — moderate confidence |
| `CALLERS=HC` only | May be a real variant rescued by joint genotyping — do not treat as low confidence |
| `CALLERS=FB` or `CALLERS=DV` only | Single independent caller — treat with caution, requires validation |

---

### 3. QUAL scores are not comparable between callers

QUAL scores from FreeBayes, HaplotypeCaller, and DeepVariant are produced by different statistical models and are not on the same scale. The merged QUAL column reflects whichever caller had the highest priority for that variant. Do not use QUAL to compare confidence across rows from different callers. Use the `CALLERS` concordance count and `FILTER=PASS` status as the primary quality signal instead.

---

### 4. Caller-specific FORMAT fields

FORMAT fields differ between callers (e.g. DeepVariant uses `VAF`; HaplotypeCaller uses `AD`, `F1R2`, `F2R1`). The merge keeps the complete FORMAT and sample column from the highest-priority caller that found the variant. Fields unique to lower-priority callers are not carried over, but the `CALLERS` tag records which other callers also found the variant.

---

### 5. Alternative tools for stricter ensemble calling

For stricter multi-caller merging consider:
- [`bcftools merge`](https://samtools.github.io/bcftools/bcftools.html#merge) — standard VCF merge with full header reconciliation
- [GLnexus](https://github.com/dnanexus-rnd/GLnexus) — joint genotyping across callers
- [GATK CombineVariants](https://gatk.broadinstitute.org/hc/en-us/articles/360037053272)

---

## License

MIT
