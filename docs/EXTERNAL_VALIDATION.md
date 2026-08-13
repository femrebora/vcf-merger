# External validation cohorts

Historical clinical cases may become an **independent** external validation cohort:

```text
historical case
  → completeness review
  → phenotype confirmation
  → variant confirmation
  → independent expert review
  → GOLD_STANDARD
  → EXTERNAL_VALIDATION dataset role
```

Entities `ExternalCohort` / `ExternalCohortVersion` exist for tracking. Membership for ML consumption still flows through versioned datasets with role `EXTERNAL_VALIDATION`.

**Invariant:** cases in EXTERNAL_VALIDATION must not later appear in TRAINING without an explicit governance operation that passes leakage checks (default: blocked).
