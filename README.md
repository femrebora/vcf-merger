# pHelix Vault

Secure genomic data **governance and storage** for the pHelix clinical genomics ecosystem.

```text
UPLOAD  ≠  PERMISSION TO TRAIN
CLINICAL STORAGE  ≠  ML TRAINING DATASET
```

Vault is the secure boundary around sensitive genomic and clinical data. It is **not** the analysis pipeline, **not** the interpretation engine, and **not** the ML training repository.

## Ecosystem position

```text
phelix-pipeline → phelix-vault → phelix-analytics → phelix-portal
                      │
                      └─ approved DatasetVersion → phelix-ml
```

## What Vault owns

- Tenant isolation
- Pseudonymous patient / identity mapping (isolated)
- Genomic asset registry (FASTQ/BAM/VCF/… metadata)
- Object storage abstraction (local + S3-compatible)
- Upload sessions + authorized short-lived download
- Processing purposes & governance metadata
- Explicit ML eligibility + curation status
- Versioned datasets with TRAINING / EXTERNAL_VALIDATION separation
- Audit events (hash-chained, append-oriented)
- Retention policy evaluation

## Quick start (local)

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
cp .env.example .env

# SQLite + local disk (default-friendly for smoke tests)
export DATABASE_URL=sqlite+pysqlite:////tmp/phelix-vault.db
export LOCAL_STORAGE_ROOT=/tmp/phelix-vault-objects
export JWT_SECRET=CHANGE_ME_DEV_ONLY_NOT_FOR_PRODUCTION

uvicorn phelix_vault.main:app --reload --port 8000
```

Or with Docker Compose (Postgres + API + MinIO):

```bash
docker compose up --build
```

OpenAPI docs: `http://localhost:8000/docs`

## Definition-of-done demo (API)

1. Mint a bootstrap admin token: `POST /api/v1/auth/dev-token`
2. Create tenant → patient → sample → case
3. `POST /api/v1/assets/upload-sessions` then complete with synthetic VCF
4. Download only with `ASSET_DOWNLOAD`
5. Set curation `TECHNICALLY_VALID`; ML remains `NOT_EVALUATED` until explicit approval
6. Approve ML eligibility; create TRAINING and EXTERNAL_VALIDATION dataset versions
7. Confirm patient overlap across conflicting roles is rejected (`dataset_leakage`)
8. Export manifest; inspect `/api/v1/audit/events`

Synthetic fixtures live under `tests/fixtures/synthetic/` only.

## Tests

```bash
pytest -q
```

## Documentation

| Doc | Topic |
|-----|--------|
| [ARCHITECTURE_ASSESSMENT.md](docs/ARCHITECTURE_ASSESSMENT.md) | Ecosystem inspection & decisions |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | System architecture |
| [SECURITY_MODEL.md](docs/SECURITY_MODEL.md) | AuthN/Z, encryption, logging |
| [THREAT_MODEL.md](docs/THREAT_MODEL.md) | Practical threats & mitigations |
| [DATA_GOVERNANCE.md](docs/DATA_GOVERNANCE.md) | Purposes, authorizations, KVKK engineering principles |
| [DATA_LIFECYCLE.md](docs/DATA_LIFECYCLE.md) | Zones, retention, deletion |
| [STORAGE_ARCHITECTURE.md](docs/STORAGE_ARCHITECTURE.md) | Object storage & keys |
| [TENANT_ISOLATION.md](docs/TENANT_ISOLATION.md) | Isolation strategy |
| [ML_DATASET_GOVERNANCE.md](docs/ML_DATASET_GOVERNANCE.md) | Eligibility & leakage controls |
| [EXTERNAL_VALIDATION.md](docs/EXTERNAL_VALIDATION.md) | Independent cohorts |
| [API.md](docs/API.md) | API groups |
| [DEPLOYMENT.md](docs/DEPLOYMENT.md) | Deploy notes |
| [DEVELOPMENT.md](docs/DEVELOPMENT.md) | Local development |
| [IMPLEMENTATION_REVIEW.md](docs/IMPLEMENTATION_REVIEW.md) | What shipped / deferred |
| [adr/](docs/adr/) | Architecture decision records |

## Important disclaimer

Technical controls support an organization’s privacy and security obligations but **do not by themselves establish legal or regulatory compliance** (including KVKK). Processing purposes, legal bases, notices, contracts, controller/processor roles, retention requirements, and international transfers require organization-specific legal/governance review.

## License

MIT
