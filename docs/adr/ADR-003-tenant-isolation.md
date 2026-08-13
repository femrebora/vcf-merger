# ADR-003: Tenant isolation strategy

## Status

Accepted

## Context

Multi-laboratory tenancy is a hard security requirement.

## Decision

Enforce tenant checks in application services for every sensitive operation; add PostgreSQL RLS later with dedicated tests.

## Consequences

Predictable fail-closed behavior now; defense-in-depth path documented.
