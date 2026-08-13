# Development

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
cp .env.example .env
pytest -q
```

## Layout

```text
src/phelix_vault/
  api/ application/ domain/ infrastructure/ config/
migrations/
tests/{unit,integration,security,fixtures/synthetic}
docs/
```

## Quality commands

```bash
ruff check src tests
ruff format src tests
mypy src/phelix_vault
pytest -q
```

## Synthetic data only

Use `tests/fixtures/synthetic/`. Never commit real patient genomic files.

## Bootstrap flow

1. `POST /api/v1/auth/dev-token` with `SYSTEM_ADMIN` and placeholder `tenant_id`
2. `POST /api/v1/tenants`
3. Mint a new token with the real `tenant_id` and `TENANT_ADMIN`
4. Exercise patient → sample → asset → case → dataset flows
