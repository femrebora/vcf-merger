# Storage architecture

## Principle

Large genomic objects are **never** stored in PostgreSQL. The database holds registry metadata; bytes live behind `ObjectStorage`.

## Interface

```text
put / get / delete / exists / checksum_sha256 / signed_download / signed_upload
```

Implementations:

- `LocalObjectStorage` — development & tests
- `S3ObjectStorage` — MinIO, AWS S3, private S3-compatible endpoints

## Object keys

```text
tenants/{tenant_id}/assets/{asset_id}/{asset_type}
```

Opaque. No patient names, national IDs, or clinical labels.

## Integrity

On finalize: size + SHA-256 computed/verified and stored. Download authorization issues short-lived URLs only.

## Encryption metadata

Assets store KMS key id + wrapped DEK fields so key rotation does not require redesigning the schema. Application-level envelope encryption of object bytes is architecture-ready; full orchestration is deferred.

## Residency

Each storage backend exposes `residency` (`TR` / `EU` / `US` / `OTHER`). Upload sessions fail if tenant `required_residency` does not match.
