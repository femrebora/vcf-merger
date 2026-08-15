# vcf-merger

Standards-aware **VCF normalization, evidence aggregation, and harmonization** for multi-caller variant callsets.

This is a **technical** ensemble / harmonization utility. It does **not** perform ACMG/AMP pathogenicity classification, AMP/ASCO/CAP somatic tiering, or any clinical interpretation. Caller concordance, QUAL, FILTER, and allele fraction are technical signals only—not evidence of clinical significance.

## What it does

- Validates input compatibility (assembly, contig dictionaries, samples, VCF vs gVCF)
- Normalizes alleles with reference-aware `bcftools norm` (left-align, trim, multiallelic split)
- Detects callers via CLI override → VCF header → filename fallback
- Aggregates **per-caller evidence** for each canonical small variant (does not discard lower-priority callers)
- Emits a compact harmonized VCF plus evidence and provenance sidecars
- Supports germline and somatic analysis modes with explicit semantics

## What it does not do

- ACMG/AMP Pathogenic / Likely Pathogenic / VUS / Likely Benign / Benign classification
- Treat caller count or concordance as pathogenicity evidence
- Silently merge GRCh37 and GRCh38 (or incompatible contig dictionaries)
- Ensemble-merge GATK HaplotypeCaller **gVCFs** (rejected by default)
- Merge structural variants / CNVs with the small-variant algorithm
- Claim clinical validation

## Supported callers (initial)

| Caller | Notes |
|--------|--------|
| GATK HaplotypeCaller | Germline VCF (not gVCF) |
| FreeBayes | Germline |
| DeepVariant | Germline |
| GATK Mutect2 | Somatic-oriented evidence |
| Strelka / Strelka2 | Germline or somatic evidence preserved natively |

Additional callers can be added via adapter modules without changing core merge logic.

## Supported file types

| Type | Support |
|------|---------|
| `.vcf` | Yes |
| `.vcf.gz` (+ tabix/CSI) | Yes |
| `.bcf` | Read where HTSlib/pysam supports it |
| `.g.vcf` / `.g.vcf.gz` | Detected and **rejected** for ensemble merge by default |

## Install

Requires Python 3.10+ and **bcftools** on `PATH` for normalization.

```bash
pip install -e .
# optional
pip install -e ".[dev]"
```

External tools: `bcftools`, `tabix` / `bgzip` (HTSlib).

## CLI

```bash
vcf-merger inspect sample.HC.vcf.gz

vcf-merger normalize \
  --reference GRCh38.fa \
  --input sample.HC.vcf.gz \
  --output sample.HC.norm.vcf.gz

vcf-merger merge \
  --mode germline \
  --reference GRCh38.fa \
  --input sample.FB.vcf.gz \
  --input sample.HC.vcf.gz \
  --input sample.DV.vcf.gz \
  --strategy union \
  --output sample.harmonized.vcf.gz
```

Somatic example:

```bash
vcf-merger merge \
  --mode somatic \
  --reference GRCh38.fa \
  --tumor-sample TUMOR \
  --normal-sample NORMAL \
  --input mutect2.vcf.gz \
  --input strelka.vcf.gz \
  --output sample.somatic.harmonized.vcf.gz
```

### Ensemble strategies

| Strategy | Behavior |
|----------|----------|
| `union` (default) | Keep alleles seen by ≥1 caller |
| `pass-union` | Keep alleles with ≥1 PASS caller |
| `consensus` | Require `--consensus-n` distinct callers (technical consensus only) |
| `caller-specific` | Include/exclude via `--include-caller` / `--exclude-caller` |

These strategies describe **technical** ensemble membership. They are not clinical filters.

## Python API

```python
from vcf_merger import harmonize_vcfs, inspect_vcf, normalize_vcf

normalize_vcf("sample.HC.vcf.gz", "sample.HC.norm.vcf.gz", reference="GRCh38.fa")
meta = inspect_vcf("sample.HC.norm.vcf.gz")
result = harmonize_vcfs(
    ["sample.FB.vcf.gz", "sample.HC.vcf.gz", "sample.DV.vcf.gz"],
    "sample.harmonized.vcf.gz",
    reference="GRCh38.fa",
    mode="germline",
    strategy="union",
)
```

## Data model

Variants and evidence are separate:

```text
CanonicalVariant (assembly, contig, normalized POS/REF/ALT)
  └── CallerEvidence[]   # never discarded merely for caller priority
        └── SampleEvidence (GT, ploidy-tolerant, DP, AD, AF, GQ, raw FORMAT)
```

There is **no** fixed scientific priority list such as MU > FB > HC > ST > DV.

## Output

For `--output sample.harmonized.vcf.gz` the tool writes:

1. `sample.harmonized.vcf.gz` (+ `.tbi` / `.csi` when possible)
2. `sample.harmonized.evidence.jsonl` — full per-caller evidence for reconstruction
3. `sample.harmonized.provenance.json` — tool version, inputs, checksums, parameters

### Harmonized VCF INFO (`VM_*`)

| Field | Meaning |
|-------|---------|
| `VM_CALLERS` | Callers that detected the allele |
| `VM_PASS_CALLERS` | Callers with FILTER=PASS |
| `VM_CALLER_COUNT` | Distinct detecting callers |
| `VM_PASS_CALLER_COUNT` | Distinct PASS callers |
| `VM_GT_CONFLICT` | Flag when genotypes disagree |
| `VM_ORIG` | One original unnormalized representation |

Absence from another caller is **not** interpreted as homozygous reference unless reference-confidence / callability data supports that conclusion (gVCF merging is out of scope).

## Germline vs somatic

- `--mode germline` — genotype comparison, allele balance, depth, GQ, phasing-friendly fields; no ACMG classification.
- `--mode somatic` — tumor/normal or tumor-only roles; preserves Mutect2/Strelka metrics natively rather than forcing a universal QUAL.

Somatic CLI:

```bash
# Matched tumor/normal
vcf-merger merge --mode somatic --reference GRCh38.fa \
  --tumor-sample TUMOR --normal-sample NORMAL \
  -i mutect2.vcf.gz -i strelka.vcf.gz -o sample.somatic.harmonized.vcf.gz

# Tumor-only
vcf-merger merge --mode somatic --reference GRCh38.fa \
  --tumor-sample TUMOR --tumor-only \
  -i mutect2.vcf.gz -o sample.tumor_only.harmonized.vcf.gz
```

Somatic INFO extras (technical only):

| Field | Meaning |
|-------|---------|
| `VM_SOMATIC_PAIR` | `tumor` or `tumor:normal` roles used |
| `VM_TUMOR_AF` | Representative tumor allele fraction |
| `VM_NORMAL_AF` | Representative normal allele fraction when available |

No AMP/ASCO/CAP clinical tier classification is performed.

Do not mix germline and somatic semantics in one run.

## Normalization

Normalization is a first-class stage (default on for `merge`). It uses `bcftools norm` for REF checks, trimming, left alignment, and multiallelic decomposition. Original representations are retained in evidence (`VM_ORIG` / sidecar).

**Tradeoff:** pysam is used for typed VCF I/O; `bcftools norm` is preferred over hand-rolled allele edits so Number=A/R/G fields are not corrupted.

## gVCF limitations

HaplotypeCaller gVCFs (`<NON_REF>`, `END` reference blocks) are intermediate reference-confidence files. Direct ensemble merging is rejected with an actionable error. Genotype gVCFs first, then merge variant-level VCFs. Extension hooks exist for future explicit gVCF support; gVCFs are never silently flattened.

## WES / WGS suitability

The harmonizer:

1. Normalizes each caller to a temporary bgzipped+indexed VCF (when `--reference` is set)
2. Streams a **multi-way merge** ordered by the reference/VCF contig dictionary
3. Writes the harmonized VCF and evidence JSONL **incrementally** (no full multi-caller DataFrame)

Peak memory is intended to stay near one contig × N callers of evidence objects, not the cartesian product of entire WGS callsets. Prefer `.vcf.gz` + tabix inputs when skipping normalization.

## Structural variants / CNVs

First-class support: SNVs and small indels (MNVs after normalization policy). Symbolic alleles (`<DEL>`, `<DUP>`, `<INV>`, `<CNV>`, BNDs) are detected and kept out of the small-variant VCF path (recorded as unsupported in the evidence sidecar). SV/CNV merging is a future extension.

## Mitochondrial / ploidy

Genotype shape is taken from the source VCF (haploid, diploid, etc.). Contig order follows the reference/VCF dictionary (not a hard-coded 1–22,X,Y,MT-only list).

## Testing

```bash
pip install -e ".[dev]"
pytest
```

Integration tests require `bcftools`. Performance tests are marked `@pytest.mark.performance`.

## Legacy scripts

`Merge_All_VCFs.py` and `Merge_all_VCF_Groups.py` remain as thin deprecated wrappers around the new API. See [docs/migration.md](docs/migration.md).

## Parser choice

**pysam (`VariantFile`)** is the primary parser/writer (HTSlib-backed typing, samples, streaming, bgzip). Pandas is not used for VCF semantics.

## Limitations

- Small-variant focus only
- gVCF ensemble merging not supported
- Representative genotype is reconciled explicitly; conflicts set `VM_GT_CONFLICT`
- QUAL scores are not comparable across callers and are not unified into one scale
- No clinical validation claim is made for this software
