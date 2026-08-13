# Tenant isolation

## Strategy (MVP)

**Primary enforcement:** application/repository filtering on `tenant_id` for every sensitive read/write, plus RBAC.

**Database RLS:** supported as a future hardening flag (`ENABLE_RLS`). Not blindly enabled without dedicated PostgreSQL integration tests.

| Layer | Status |
|-------|--------|
| API principal `tenant_id` claim | Required |
| Service/repository filters | Required |
| Identical not-found across tenants | Required |
| Automated cross-tenant tests | Required |
| PostgreSQL RLS | Documented / deferred |

## Why not RLS-only

RLS is valuable defense-in-depth but easy to misconfigure with migrations, admin roles, and background jobs. Vault starts with explicit application checks that are unit/integration tested, then adds RLS once policies and test harnesses are ready.
