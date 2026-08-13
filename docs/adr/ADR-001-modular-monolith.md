# ADR-001: Modular monolith

## Status

Accepted

## Context

pHelix Vault needs strong domain boundaries without early microservice operational cost.

## Decision

Ship a single deployable service with package-level domains (identity, assets, governance, datasets, audit).

## Consequences

Clear invariants and simpler local development; physical split remains possible later if scale/isolation demands it.
