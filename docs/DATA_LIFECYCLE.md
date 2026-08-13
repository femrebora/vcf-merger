# Data lifecycle

## Logical zones

```text
RAW_CLINICAL → PROCESSED_CLINICAL → CURATED → ML_APPROVED
```

Zones are metadata in MVP and may map to separate physical buckets later.

## Asset lifecycle

```text
UPLOAD_PENDING → AVAILABLE → PENDING_DELETION → DELETED
```

Deletion removes object bytes when the storage adapter can, records an audit event, and does **not** claim guaranteed physical erasure beyond provider capabilities.

## Curation status (explicit)

```text
RAW → TECHNICALLY_VALID → … → GOLD_STANDARD
```

Non-sequential jumps require an explicit reason. Transitions are stored with actor, timestamps, and prior state.

## ML eligibility (explicit)

```text
NOT_EVALUATED → PENDING_REVIEW → APPROVED | INELIGIBLE | REVOKED
```

`REVOKED` blocks inclusion in future dataset versions. Historical approved dataset versions remain immutable for reproducibility.

## Retention

`RetentionPolicy` defines period + action (`REVIEW` / `ARCHIVE` / `DELETE`).  
`evaluate` reports due assets; **automatic deletion is not enabled by default in development**.
