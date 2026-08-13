# pHelix Vault — Architecture Assessment

**Date:** 2026-08-13  
**Workspace inspected:** `github.com/femrebora/vcf-merger` (agent run titled “pHelix genomic data vault”)  
**Branch for this work:** `cursor/phelix-vault-foundation-f2cd`

---

## 1. Existing pHelix conventions

### Finding: no sibling `phelix-*` repositories are available

Searched the authenticated GitHub account and workspace for:

- `phelix-pipeline`
- `phelix-analytics`
- `phelix-portal`
- `phelix-vault`
- `phelix-ml`

**Result:** none exist in this workspace or under the accessible GitHub account. There is therefore **no established multi-repo pHelix API/Docker/auth convention to inherit**.

### Closest available codebase (current workspace)

The attached repository is **`vcf-merger`**, a Python bioinformatics library for technical VCF harmonization. Observed conventions:

| Area | Convention in `vcf-merger` |
|------|----------------------------|
| Language | Python ≥3.10 |
| Packaging | `pyproject.toml` + setuptools (`src/` layout) |
| Testing | `pytest`, `tests/unit`, `tests/integration`, `tests/fixtures` |
| CI | GitHub Actions: checkout → setup-python 3.11 → `pip install -e ".[dev]"` → `pytest` |
| Typing / models | `dataclasses`, `Enum`, `Protocol`, `__future__` annotations |
| Logging | stdlib `logging` in CLI |
| Docker | **none** |
| Database | **none** |
| Auth | **none** |
| API | CLI + library API only (no HTTP) |
| Domain philosophy | Explicit “no clinical interpretation” boundaries; provenance sidecars; synthetic fixtures |

Related public forks under the same account (not pHelix products, but relevant tooling interest): MinIO, MISO LIMS, OpenELIS Global — suggest eventual comfort with **S3-compatible object storage** and LIMS-style clinical systems, but they are not a shared pHelix stack.

### Implication for Vault

pHelix Vault is a **greenfield service**. We align with the only local convention that exists: **Python + `src/` layout + pytest + GitHub Actions**, and introduce FastAPI/PostgreSQL/object-storage patterns appropriate for a secure data boundary service.

This branch **re-targets the working tree toward pHelix Vault**. `vcf-merger` remains in git history on `main`; it is not the Vault product and must not be merged into Vault responsibilities.

---

## 2. Technologies detected

| Layer | Detected | Proposed for Vault |
|-------|----------|--------------------|
| Language | Python 3.10+ / CI 3.11 | Python 3.11 |
| Package manager | pip / setuptools / `pyproject.toml` | Same (`pyproject.toml`) |
| HTTP API | None | **FastAPI** + Pydantic v2 (OpenAPI, typed DTOs) |
| DB | None | **PostgreSQL** + SQLAlchemy 2.0 + Alembic |
| Object storage | None (MinIO fork interest) | **Local FS** + **S3-compatible** (MinIO/AWS) |
| Auth | None | **JWT bearer** with OIDC-ready claims; dev tokens for local |
| Containers | None | Docker + Compose (`vault-api`, `postgres`, `minio`) |
| Lint/types | Minimal | `ruff`, `mypy`, `pytest` |
| CI | Simple pytest job | Extended: lint, typecheck, unit/integration/security tests, migration check |

---

## 3. Assumptions

1. **Greenfield ecosystem:** Vault defines initial HTTP, authz, and storage conventions that future `phelix-*` services should follow unless a later shared platform repo supersedes them.
2. **Modular monolith:** One deployable API process with clear domain packages — not a microservice mesh.
3. **PostgreSQL** is the system of record for metadata; genomic blobs live in object storage.
4. **Authentication** is externalizable (JWT / future OIDC). Vault owns **authorization** (RBAC).
5. **Development auth** may use signed HS256 JWTs from an env secret — not a homegrown password DB.
6. **Tenant isolation** is enforced in repositories + policy layer; PostgreSQL RLS is designed but optional behind a flag until thoroughly tested.
7. **No real patient data** in repo; only synthetic fixtures under `tests/fixtures/synthetic/`.
8. Working in `vcf-merger` remote is a **workspace attachment mismatch**; deliverable is Vault. Prefer creating/renaming a dedicated `phelix-vault` GitHub repository before production use.

---

## 4. Integration boundaries

```
phelix-pipeline  →  issues temporary asset access; returns derived asset refs + checksums
phelix-vault     →  owns identity mapping, genomic assets, governance, datasets, audit
phelix-analytics →  receives case/sample/VCF/phenotype refs via controlled APIs
phelix-portal    →  never talks to object storage; always via Vault APIs
phelix-ml        →  consumes versioned Dataset manifests only; never raw SQL over Vault
```

Vault **must not**:

- run alignment / variant calling
- perform ACMG/AMP (or other) clinical interpretation
- train ML models
- claim KVKK/GDPR “compliance” as a software property

---

## 5. Proposed pHelix Vault architecture

### Shape

Modular monolith:

```
API (FastAPI) → Application services → Domain rules → Infrastructure adapters
```

### Domain packages

| Domain | Owns |
|--------|------|
| `tenants` | Tenant, membership, residency policy |
| `identity` | PseudonymousPatient, IdentityMapping (isolated) |
| `samples` | Sample registry |
| `assets` | GenomicAsset, upload sessions, provenance refs |
| `cases` | GenomicCase, curation status, data-quality dimensions |
| `governance` | ProcessingPurpose, ProcessingAuthorization, ML eligibility |
| `datasets` | Dataset, DatasetVersion, DatasetCase, leakage checks, ExternalCohort |
| `retention` | RetentionPolicy, deletion lifecycle |
| `audit` | Append-oriented AuditEvent records |

### Storage / crypto

- `ObjectStorage` protocol: Local + S3-compatible
- `KeyManagementProvider`: DevelopmentKeyProvider now; ExternalKMSProvider stub
- Envelope-encryption-ready metadata on assets (`key_id`, `wrapping_alg`, etc.)

### AuthZ

Central `AuthorizationService` mapping roles → permissions. Handlers call `require(permission, resource)` — no scattered role `if`s.

### Invariants (enforced in domain + tests)

1. Sensitive objects always carry `tenant_id`
2. Identity mapping not required for ordinary asset/analysis flows
3. Upload ≠ ML eligibility
4. EXTERNAL_VALIDATION cannot silently enter TRAINING
5. Downloads authorized + audited
6. No genomic blobs in git
7. No permanent public URLs
8. No genomic/PHI content in app logs
9. Dataset versions used for ML are immutable membership snapshots
10. Vault never interprets variants clinically

---

## 6. Security assumptions

| Topic | Assumption |
|-------|------------|
| Trust boundary | Vault API is the only path to object bytes for portal/analytics/ML |
| Secrets | Env / runtime injection only; `.env.example` placeholders |
| Transport | TLS terminated at reverse proxy in production; local HTTP OK for compose |
| At-rest | Object storage encryption + optional app-level DEK metadata |
| Audit | Separate from operational logs; no PHI payloads |
| Cross-tenant | Denied by default; security tests must fail closed |
| Legal | Structures record decisions; software does not assert compliance |

---

## 7. Unanswered questions

1. Will a dedicated `femrebora/phelix-vault` (or org) repository be created, or should this continue under `vcf-merger`?
2. Preferred production IdP (Keycloak, Auth0, hospital AD OIDC, etc.)?
3. Required initial data residency default (`TR` vs configurable per tenant)?
4. Will identity mappings live in the same PostgreSQL instance as asset metadata in v1, or a split “identity vault” DB immediately?
5. Existing hospital object-storage endpoint (MinIO / cloud) and network constraints?
6. Retention default periods and who may approve deletion?
7. Should portal and analytics share the same JWT audience, or use service credentials?

Until answered, Vault uses: same DB with isolated identity tables, residency metadata defaulting to tenant config (`TR` example), and JWT roles embedded in claims for development.

---

## 8. Implementation phases

### Phase A — Assessment (this document)

Complete.

### Phase B — Foundation MVP (this delivery)

- Domain models + SQLAlchemy mappings + Alembic migrations
- Tenant / patient / sample / asset / case / governance / dataset / audit / retention
- Local + S3 storage adapters
- Upload-session + finalize + authorized download
- RBAC policy layer + JWT auth
- Dataset versioning + TRAINING vs EXTERNAL_VALIDATION leakage checks
- Docker Compose (api, postgres, minio)
- Unit / integration / security tests + synthetic fixtures
- Core docs + ADRs + threat model
- CI pipeline

### Phase C — Hardening (near-term follow-on)

- PostgreSQL RLS with comprehensive tests
- External KMS integration
- VCF header sanitization pipeline (ORIGINAL / SANITIZED linkage)
- Immutable / WORM audit sink
- Rate limiting / SIEM hooks
- Full legal-hold automation

### Phase D — Ecosystem wiring

- Pipeline run handshake APIs
- Analytics limited-scope grant tokens
- Portal UX contracts
- ML dataset export tokens + reproducibility records

---

## Decision summary

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Service shape | Modular monolith | MVP clarity; strong domains without ops overhead |
| Language | Python 3.11 | Matches available ecosystem tooling |
| API | FastAPI | OpenAPI, typing, async I/O for uploads |
| DB | PostgreSQL | Preferred in mission; no conflicting standard |
| Blobs | Object storage | Never store FASTQ/BAM/VCF in Postgres |
| AuthN | JWT (OIDC-ready) | No homegrown password store |
| AuthZ | Central RBAC | Room for ABAC later |
| ML boundary | Explicit eligibility + dataset versions | Upload ≠ permission to train |
