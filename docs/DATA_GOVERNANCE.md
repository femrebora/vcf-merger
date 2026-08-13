# Data governance

Vault records **governance metadata**. It does not make legal determinations.

## Processing purposes

Explicit purposes include clinical diagnosis, reanalysis, QC, validation, research, ML training, and ML external validation. Purposes must be set deliberately — never inferred from upload alone.

## ProcessingAuthorization

Stores configurable fields such as legal basis code, authorization source, effective dates, restrictions, controller/processor, international transfer flag, and ML/research/reanalysis allowances.

Software must not claim “KVKK compliant” merely because these fields exist.

## KVKK-oriented engineering principles

Vault implements engineering support for:

- purpose limitation
- data minimization (opaque IDs; identity isolation)
- authorization
- auditability
- retention evaluation
- pseudonymization
- special-category data protection posture (treat genomics as highly sensitive)
- data residency metadata
- international-transfer controls (flags + residency checks)
- ML-use restrictions

**Disclaimer:** Technical controls support an organization’s privacy and security obligations but do not by themselves establish legal or regulatory compliance. Processing purposes, legal bases, notices, contracts, controller/processor roles, retention requirements, and international transfers require organization-specific legal/governance review.
