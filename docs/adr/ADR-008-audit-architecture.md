# ADR-008: Audit architecture

## Status

Accepted

## Context

Sensitive actions require accountability without logging PHI payloads.

## Decision

Append-oriented `AuditEvent` table with hash chaining; separate from operational logs; ban genomic/identity payloads in metadata.

## Consequences

Tamper-evident readiness now; WORM/SIEM sinks later.
