# ADR-002: Object storage vs database blobs

## Status

Accepted

## Context

FASTQ/BAM/VCF objects are multi-gigabyte.

## Decision

Store bytes in object storage; metadata and checksums in PostgreSQL.

## Consequences

Streaming uploads/downloads; portable S3 interface; no Postgres bloat.
