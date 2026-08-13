# Deployment

## Containers

`docker compose up --build` starts:

- `vault-api` (FastAPI / Uvicorn)
- `postgres`
- `minio` (+ bucket init)

Do not bake production secrets into images. Inject via environment or a secret manager.

## Required production settings

- `ENVIRONMENT=production` (disables dev token route)
- Strong `JWT_SECRET` or preferably external OIDC validation
- `DATABASE_URL` → PostgreSQL
- `STORAGE_BACKEND=s3` with scoped credentials
- Matching `STORAGE_RESIDENCY` and tenant policies
- TLS at the ingress / reverse proxy
- Restrict `CORS_ALLOWED_ORIGINS`

## Migrations

```bash
alembic upgrade head
```

MVP also calls `create_all` on startup for local ergonomics; production should prefer Alembic-only schema management.

## Health

- Liveness: `GET /health`
- Readiness: `GET /ready`
