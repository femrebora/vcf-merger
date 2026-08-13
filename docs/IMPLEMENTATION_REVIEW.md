# Implementation review

## What was implemented

- Modular FastAPI service with domain-separated SQLAlchemy models
- Tenant, pseudonymous patient, identity mapping, sample, genomic case
- Genomic asset registry + upload session + checksum finalize + authorized download
- Local and S3-compatible object storage adapters
- Development KMS + external KMS stub
- Central RBAC + JWT auth (dev token route)
- Processing purposes, ProcessingAuthorization, retention policies
- Explicit curation + ML eligibility transitions with history
- Dataset / DatasetVersion / manifest export
- TRAINING vs EXTERNAL_VALIDATION (and TEST) patient-leakage guards
- Hash-chained audit events
- Alembic initial migration, Docker Compose, CI, synthetic fixtures
- Architecture, security, governance, ADR, and threat-model documentation

## Intentionally deferred

- External enterprise KMS orchestration
- PostgreSQL RLS enablement with full test matrix
- Legal hold automation / consent UI
- Advanced ABAC policy engine
- Physical multi-region routing
- Immutable WORM audit backend / SIEM / DLP
- Full VCF header sanitization worker (ORIGINAL/SANITIZED linkage modeled)
- Direct browser→S3 multipart for multi-hundred-GB production uploads
- Family / variant-level leakage analytics beyond patient overlap

These are marked as stubs or docs — not faked as complete security features.

## Security controls implemented

- Tenant checks on sensitive reads/writes
- Permission checks via centralized RBAC
- Opaque IDs and safe object keys
- Short-lived download URLs
- Audit for upload/download/eligibility/dataset events
- Safe logging redaction helpers
- Secret-free `.env.example` + gitignore for genomic artifacts
- Residency mismatch rejection at upload session creation

## Limitations / known risks

- Workspace was attached to `femrebora/vcf-merger`; no sibling `phelix-*` repos existed. Stack chosen as greenfield Python aligned with available tooling.
- Application-level tenant isolation is primary; RLS not yet enforced.
- Dev JWT secret must never be used in production.
- Local finalize path streams via spooled temp file; production should prefer direct-to-object-storage multipart.
- Identity mappings currently share the same DB instance (table isolation only).

## Recommended next steps

1. Create/rename a dedicated `phelix-vault` GitHub repository
2. Wire real OIDC provider; remove shared HS256 secret
3. Enable Postgres RLS with integration tests
4. Implement VCF metadata sanitization pipeline
5. Add pipeline handshake APIs for `phelix-pipeline`
6. Add scoped export tokens for `phelix-ml` and analytics

## Integration points

| System | Contract |
|--------|----------|
| phelix-pipeline | Temporary input access + register derived assets/checksums/versions |
| phelix-analytics | Case/sample/VCF/phenotype references via controlled APIs |
| phelix-portal | All access via Vault API; never direct object storage |
| phelix-ml | Consume approved DatasetVersion manifests only |
