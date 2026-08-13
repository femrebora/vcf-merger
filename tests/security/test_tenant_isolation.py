"""Security tests: cross-tenant and authorization boundaries."""

from __future__ import annotations

import hashlib

from fastapi.testclient import TestClient

from tests.helpers import auth_header, bootstrap_tenant


def _seed_asset(client: TestClient, headers: dict[str, str]) -> str:
    patient = client.post("/api/v1/patients", headers=headers, json={}).json()
    sample = client.post(
        "/api/v1/samples",
        headers=headers,
        json={"patient_id": patient["id"]},
    ).json()
    content = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    digest = hashlib.sha256(content).hexdigest()
    session = client.post(
        "/api/v1/assets/upload-sessions",
        headers=headers,
        json={
            "patient_id": patient["id"],
            "sample_id": sample["id"],
            "asset_type": "VCF",
            "content_type": "text/plain",
        },
    ).json()
    asset = client.post(
        f"/api/v1/assets/upload-sessions/{session['session_id']}/complete",
        headers=headers,
        files={"file": ("x.vcf", content, "text/plain")},
        data={"expected_sha256": digest},
    ).json()
    return asset["id"]


def test_cross_tenant_asset_access_denied(client: TestClient) -> None:
    _, headers_a = bootstrap_tenant(client, "Tenant A")
    _, headers_b = bootstrap_tenant(client, "Tenant B")
    asset_id = _seed_asset(client, headers_a)

    meta = client.get(f"/api/v1/assets/{asset_id}", headers=headers_b)
    assert meta.status_code == 404

    download = client.post(f"/api/v1/assets/{asset_id}/download", headers=headers_b)
    assert download.status_code == 404


def test_ml_researcher_cannot_read_identity(client: TestClient) -> None:
    tenant_id, admin = bootstrap_tenant(client, "Tenant Identity")
    patient = client.post(
        "/api/v1/patients",
        headers=admin,
        json={"external_subject_id": "HOSP-SYNTH-999"},
    ).json()

    ml = auth_header(client, tenant_id=tenant_id, roles=["ML_RESEARCHER"], subject="ml-1")
    resp = client.get(f"/api/v1/patients/{patient['id']}/identity-mapping", headers=ml)
    assert resp.status_code == 403


def test_auditor_cannot_download(client: TestClient) -> None:
    tenant_id, admin = bootstrap_tenant(client, "Tenant Auditor")
    asset_id = _seed_asset(client, admin)
    auditor = auth_header(client, tenant_id=tenant_id, roles=["AUDITOR"], subject="aud-1")
    resp = client.post(f"/api/v1/assets/{asset_id}/download", headers=auditor)
    assert resp.status_code == 403


def test_expired_upload_session_rejected(client: TestClient, monkeypatch) -> None:
    from datetime import timedelta

    from phelix_vault.domain import common as common_mod

    _, headers = bootstrap_tenant(client, "Tenant Expiry")
    patient = client.post("/api/v1/patients", headers=headers, json={}).json()
    sample = client.post(
        "/api/v1/samples",
        headers=headers,
        json={"patient_id": patient["id"]},
    ).json()
    session = client.post(
        "/api/v1/assets/upload-sessions",
        headers=headers,
        json={
            "patient_id": patient["id"],
            "sample_id": sample["id"],
            "asset_type": "VCF",
            "content_type": "text/plain",
        },
    ).json()

    real_utcnow = common_mod.utcnow

    def future():
        return real_utcnow() + timedelta(hours=5)

    monkeypatch.setattr("phelix_vault.application.services.assets.utcnow", future)
    resp = client.post(
        f"/api/v1/assets/upload-sessions/{session['session_id']}/complete",
        headers=headers,
        files={"file": ("x.vcf", b"x", "text/plain")},
    )
    assert resp.status_code == 422


def test_unsafe_object_key_rejected() -> None:
    import pytest

    from phelix_vault.domain.errors import ValidationError
    from phelix_vault.infrastructure.storage.base import validate_object_key

    with pytest.raises(ValidationError):
        validate_object_key("../../secret")
