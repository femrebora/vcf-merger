# Threat model

Practical threats for pHelix Vault. This is an engineering document, not a certification.

| Threat | Attack path | Mitigation | Remaining risk |
|--------|-------------|------------|----------------|
| Cross-tenant access | Guess/leak asset IDs; query without tenant filter | Repository tenant checks; identical 404; security tests | Bug in a new query path |
| Stolen credentials | Token theft from workstation | Short JWT TTL; TLS; future OIDC + refresh | Endpoint malware |
| Developer access to prod data | Shared prod credentials in laptop env | Separate envs; no prod secrets in repo; least privilege | Insider with prod role |
| Object-storage credential leakage | Committed `.env` / CI log | `.gitignore`; placeholders only; IAM scoping | Misconfigured cloud IAM |
| SQL injection | String-built SQL | SQLAlchemy bound parameters | Raw SQL additions |
| Path traversal | `../` object keys | `validate_object_key` | Alternate encodings |
| Malicious uploads | Huge/malformed files | Size limits; content-type recorded; no exec of uploads | Zip bombs / resource abuse |
| Logs leaking genomics | `logger.info(vcf)` | Safe logging utils; code review; tests for redaction | Verbose third-party libs |
| Backup leakage | Unencrypted DB/object backups | Encrypt backups; access control (ops) | Backup vendor compromise |
| Excessive staff privileges | Over-broad TENANT_ADMIN | Role matrix; AUDITOR without download | Social pressure to elevate |
| ML dataset misuse | Train on clinical dump | Explicit eligibility + versioned manifests | Process bypass outside Vault |
| Accidental external transfer | Wrong residency backend | Tenant residency vs storage residency check | Manual ops copy |
| Dependency compromise | Malicious PyPI package | Pin/CI; future SCA scanning | Supply chain |
| Insider threat | Privileged user exfiltrates | Audit downloads; future DLP/SIEM | Determined insider |

## Abuse cases covered by automated tests

- Cross-tenant asset meta/download denied
- ML researcher blocked from identity mapping
- Auditor blocked from download
- Expired upload session rejected
- TRAINING vs EXTERNAL_VALIDATION patient overlap rejected
