# ADR-004: Pseudonymous identifiers

## Status

Accepted

## Context

Human-readable filenames and sequential IDs leak identity.

## Decision

Opaque prefixed ULIDs (`PHX-PAT`, `PHX-SMP`, `PHX-AST`, …) with isolated `IdentityMapping`.

## Consequences

Analysis can proceed without identity; mapping access is a privileged, audited action.
