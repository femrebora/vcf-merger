"""Shared test helpers (not fixtures)."""

from __future__ import annotations

from fastapi.testclient import TestClient


def auth_header(
    client: TestClient, *, tenant_id: str, roles: list[str], subject: str = "user-1"
) -> dict[str, str]:
    resp = client.post(
        "/api/v1/auth/dev-token",
        json={"subject": subject, "tenant_id": tenant_id, "roles": roles},
    )
    assert resp.status_code == 200, resp.text
    token = resp.json()["access_token"]
    return {"Authorization": f"Bearer {token}"}


def bootstrap_tenant(client: TestClient, name: str = "Lab A") -> tuple[str, dict[str, str]]:
    bootstrap = auth_header(
        client, tenant_id="bootstrap", roles=["SYSTEM_ADMIN"], subject="bootstrap"
    )
    resp = client.post(
        "/api/v1/tenants",
        headers=bootstrap,
        json={"name": name, "required_residency": "TR"},
    )
    assert resp.status_code == 200, resp.text
    tenant_id = resp.json()["id"]
    headers = auth_header(client, tenant_id=tenant_id, roles=["TENANT_ADMIN"], subject="admin-a")
    return tenant_id, headers
