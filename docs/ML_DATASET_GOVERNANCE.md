# ML dataset governance

## Core rule

```text
patient genomic data
   → secure clinical storage
   → curation
   → governance review
   → explicit ML eligibility
   → versioned dataset
   → phelix-ml
```

Upload alone never grants training rights.

## Dataset roles

`TRAINING` | `VALIDATION` | `TEST` | `EXTERNAL_VALIDATION` | `BENCHMARK` | `NONE`

## Leakage safeguards

Approved dataset versions with conflicting roles may not share patients:

- TRAINING ↔ EXTERNAL_VALIDATION
- TRAINING ↔ TEST
- VALIDATION ↔ EXTERNAL_VALIDATION

Checks run at version creation and again at approval. Violations raise `dataset_leakage`.

## Manifests

ML consumes an exported immutable-ish manifest for a specific `DatasetVersion` — never `SELECT * FROM patients`.

## Model linkage

`ModelDatasetLink` records that DatasetVersion X was used by ModelVersion Y without importing ML libraries.
