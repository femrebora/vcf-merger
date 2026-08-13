# API

Base prefix: `/api/v1`  
OpenAPI: `/docs` · `/openapi.json`

## Groups

| Group | Purpose |
|-------|---------|
| `/auth/dev-token` | Development JWT minting (not in production) |
| `/tenants` | Tenant registry |
| `/patients` | Pseudonymous patients; identity mapping (restricted) |
| `/samples` | Sample registry |
| `/cases` | Cases, curation, ML eligibility, eligibility report |
| `/assets` | Upload sessions, metadata, authorized download |
| `/datasets` | Datasets, versions, approve, manifest, overlap |
| `/governance/authorizations` | ProcessingAuthorization records |
| `/retention` | Policies + evaluation |
| `/audit/events` | Audit trail |

Health (unprefixed): `/health`, `/ready`

## Conventions

- Typed Pydantic request/response DTOs
- Domain errors → `{error, message, request_id}`
- Bearer JWT required for sensitive routes
- Correlation via `X-Request-ID`
