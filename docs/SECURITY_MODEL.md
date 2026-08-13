# Security model

## Authentication

- JWT bearer tokens (HS256 in development)
- Claims: `sub`, `tenant_id`, `roles`, `aud`, `iss`
- OIDC/OAuth2 integration is the intended production path; Vault verifies identity and enforces its own authorization
- `POST /api/v1/auth/dev-token` is **disabled when `ENVIRONMENT=production`**

## Authorization (RBAC)

Central `AuthorizationService` maps roles → permissions. Controllers must not scatter ad-hoc role checks.

| Role | Notable restrictions |
|------|----------------------|
| ML_RESEARCHER | No identity mapping access |
| LAB_TECHNICIAN | Cannot approve ML datasets |
| AUDITOR | Audit read; cannot download genomic assets |
| DATA_CURATOR | Curate cases; identity decrypt not granted by default |

Room remains for future ABAC (purpose, residency, zone).

## Encryption

- TLS in transit (terminate at reverse proxy in production)
- Object storage encryption at rest (provider-dependent)
- Envelope-encryption-ready asset metadata: `encryption_key_id`, `encrypted_data_key`, `encryption_algorithm`
- `DevelopmentKeyProvider` for local only; `ExternalKMSProvider` stub for future KMS

Never commit real keys. Use environment / secret injection.

## Logging

- Operational logs vs security audit logs are separate
- `safe_extra` / redaction helpers block sensitive keys
- Audit events are hash-chained (`prev_hash` / `event_hash`) for tamper-evidence readiness
- Audit metadata must not contain VCF/FASTQ contents, names, or national IDs

## HTTP hardening

- Restricted CORS origins
- Security headers (`nosniff`, `DENY` frame, `no-store`)
- Upload size ceiling
- Object-key path traversal rejection
- Short-lived signed download URLs only
