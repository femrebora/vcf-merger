# Architecture

pHelix Vault is a **modular monolith**: one deployable API with strong domain packages.

```text
                    pHelix Vault
            ┌─────────────────────────┐
            │       API Layer         │
            └────────────┬────────────┘
                         │
            ┌────────────▼────────────┐
            │   Application Services  │
            └────────────┬────────────┘
   ┌─────────────────────┼─────────────────────┐
   ▼                     ▼                     ▼
 Identity            Genomic Assets         Governance
 pseudonyms          object metadata         purposes
 mapping             checksums               ML eligibility
                     provenance              datasets
                         │
                         ▼
                   Audit / Security
```

## Layers

| Layer | Responsibility |
|-------|----------------|
| `api/` | HTTP routes, DTOs, auth dependency injection |
| `application/services/` | Use cases & invariants |
| `domain/` | Enums, errors, shared primitives |
| `infrastructure/` | DB, object storage, KMS, JWT, logging |

## Architectural invariants

1. Every sensitive object has `tenant_id`
2. Identity mapping is not required for ordinary genomic analysis
3. Clinical upload ≠ ML eligibility
4. EXTERNAL_VALIDATION cannot silently enter TRAINING
5. Downloads are authorized and audited
6. Raw genomic files never live in Git
7. No permanent public genomic URLs
8. No genomic/PHI contents in ordinary logs
9. Dataset versions used for ML are reproducible membership snapshots
10. Vault never performs clinical interpretation

## Identity vs genomic data

- `PseudonymousPatient` / `Sample` / `GenomicAsset` use opaque IDs (`PHX-PAT-…`)
- `IdentityMapping` is a separate table and permission (`IDENTITY_MAPPING_*`)
- Object keys are tenant/asset scoped — never patient names or national IDs

## Storage

Metadata → PostgreSQL (SQLite supported for tests).  
Blobs → `ObjectStorage` (local or S3-compatible).

## Future physical isolation

Logical data zones (`RAW_CLINICAL` → `ML_APPROVED`) are metadata today and can map to separate buckets/credentials later (IDENTITY / GENOMIC / ML vaults).
