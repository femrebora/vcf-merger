# ADR-005: Dataset versioning

## Status

Accepted

## Context

ML reproducibility requires stable membership.

## Decision

Immutable approved `DatasetVersion` snapshots exported as manifests; changes create new versions.

## Consequences

Training/evaluation always cites dataset id + version.
