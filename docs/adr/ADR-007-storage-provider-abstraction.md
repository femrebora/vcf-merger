# ADR-007: Storage-provider abstraction

## Status

Accepted

## Context

Cloud/location choices must remain configurable.

## Decision

`ObjectStorage` protocol with Local and S3-compatible adapters; residency metadata on the backend.

## Consequences

MinIO/AWS/private S3 without rewriting domain logic.
