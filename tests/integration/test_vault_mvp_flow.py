"""End-to-end MVP flow against the definition of done."""

from __future__ import annotations

import hashlib
from pathlib import Path

from fastapi.testclient import TestClient

from tests.helpers import bootstrap_tenant


def test_mvp_definition_of_done(client: TestClient, synthetic_dir: Path) -> None:
    tenant_id, headers = bootstrap_tenant(client, "Clinical Lab Alpha")

    # Pseudonymous patient (no human-readable storage names)
    patient = client.post("/api/v1/patients", headers=headers, json={}).json()
    assert patient["id"].startswith("PHX-PAT-")

    sample = client.post(
        "/api/v1/samples",
        headers=headers,
        json={"patient_id": patient["id"], "sample_type": "WES"},
    ).json()
    assert sample["id"].startswith("PHX-SMP-")

    case = client.post(
        "/api/v1/cases",
        headers=headers,
        json={
            "patient_id": patient["id"],
            "sample_id": sample["id"],
            "processing_purposes": ["CLINICAL_DIAGNOSIS"],
        },
    ).json()
    assert case["curation_status"] == "RAW"
    assert case["ml_eligibility"] == "NOT_EVALUATED"

    vcf_path = synthetic_dir / "sample.vcf"
    content = vcf_path.read_bytes()
    digest = hashlib.sha256(content).hexdigest()

    upload = client.post(
        "/api/v1/assets/upload-sessions",
        headers=headers,
        json={
            "patient_id": patient["id"],
            "sample_id": sample["id"],
            "case_id": case["id"],
            "asset_type": "VCF",
            "content_type": "text/plain",
            "genome_build": "GRCh38",
            "pipeline_name": "synthetic",
            "pipeline_version": "0.0.1",
        },
    )
    assert upload.status_code == 200, upload.text
    session = upload.json()
    assert session["asset_id"].startswith("PHX-AST-")

    complete = client.post(
        f"/api/v1/assets/upload-sessions/{session['session_id']}/complete",
        headers=headers,
        files={"file": ("sample.vcf", content, "text/plain")},
        data={"expected_sha256": digest},
    )
    assert complete.status_code == 200, complete.text
    asset = complete.json()
    assert asset["sha256"] == digest
    assert asset["status"] == "AVAILABLE"
    assert asset["size_bytes"] == len(content)

    meta = client.get(f"/api/v1/assets/{asset['id']}", headers=headers)
    assert meta.status_code == 200

    download = client.post(f"/api/v1/assets/{asset['id']}/download", headers=headers)
    assert download.status_code == 200
    assert "download_url" in download.json()

    # Curation → TECHNICALLY_VALID; ML still NOT_EVALUATED until explicit approval
    curated = client.post(
        f"/api/v1/cases/{case['id']}/curation-status",
        headers=headers,
        json={"status": "TECHNICALLY_VALID"},
    ).json()
    assert curated["curation_status"] == "TECHNICALLY_VALID"
    assert curated["ml_eligibility"] == "NOT_EVALUATED"

    approved = client.post(
        f"/api/v1/cases/{case['id']}/ml-eligibility",
        headers=headers,
        json={"status": "APPROVED", "reason": "synthetic governance approval"},
    ).json()
    assert approved["ml_eligibility"] == "APPROVED"

    # TRAINING dataset
    ds_train = client.post(
        "/api/v1/datasets",
        headers=headers,
        json={"name": "Train Set", "dataset_role": "TRAINING"},
    ).json()
    ver_train = client.post(
        f"/api/v1/datasets/{ds_train['id']}/versions",
        headers=headers,
        json={
            "version": "1.0.0",
            "purpose": "ML_TRAINING",
            "case_ids": [case["id"]],
            "eligibility_criteria": {"synthetic": True},
        },
    ).json()
    approved_train = client.post(
        f"/api/v1/datasets/versions/{ver_train['id']}/approve",
        headers=headers,
    )
    assert approved_train.status_code == 200, approved_train.text

    # Separate EXTERNAL_VALIDATION dataset with a different patient
    patient_b = client.post("/api/v1/patients", headers=headers, json={}).json()
    sample_b = client.post(
        "/api/v1/samples",
        headers=headers,
        json={"patient_id": patient_b["id"], "sample_type": "WES"},
    ).json()
    case_b = client.post(
        "/api/v1/cases",
        headers=headers,
        json={"patient_id": patient_b["id"], "sample_id": sample_b["id"]},
    ).json()
    client.post(
        f"/api/v1/cases/{case_b['id']}/ml-eligibility",
        headers=headers,
        json={"status": "APPROVED", "reason": "external cohort"},
    )

    ds_ext = client.post(
        "/api/v1/datasets",
        headers=headers,
        json={"name": "External Val", "dataset_role": "EXTERNAL_VALIDATION"},
    ).json()
    ver_ext = client.post(
        f"/api/v1/datasets/{ds_ext['id']}/versions",
        headers=headers,
        json={
            "version": "1.0.0",
            "purpose": "ML_EXTERNAL_VALIDATION",
            "case_ids": [case_b["id"]],
        },
    )
    assert ver_ext.status_code == 200, ver_ext.text
    client.post(f"/api/v1/datasets/versions/{ver_ext.json()['id']}/approve", headers=headers)

    # Same patient cannot enter conflicting EXTERNAL_VALIDATION after TRAINING approval
    leak = client.post(
        f"/api/v1/datasets/{ds_ext['id']}/versions",
        headers=headers,
        json={
            "version": "1.1.0",
            "purpose": "ML_EXTERNAL_VALIDATION",
            "case_ids": [case["id"]],
        },
    )
    assert leak.status_code == 409, leak.text
    assert leak.json()["error"] == "dataset_leakage"

    manifest = client.get(
        f"/api/v1/datasets/versions/{ver_train['id']}/manifest",
        headers=headers,
    )
    assert manifest.status_code == 200
    body = manifest.json()
    assert body["dataset_role"] == "TRAINING"
    assert body["version"] == "1.0.0"
    assert body["cases"][0]["case_id"] == case["id"]
    assert body["cases"][0]["vcf_asset_id"] == asset["id"]

    audit = client.get("/api/v1/audit/events", headers=headers)
    assert audit.status_code == 200
    actions = {e["action"] for e in audit.json()}
    assert "ASSET_UPLOADED" in actions
    assert "ASSET_DOWNLOADED" in actions
    assert "ML_ELIGIBILITY_CHANGED" in actions
    assert "DATASET_APPROVED" in actions


def test_upload_not_ml_permission(client: TestClient) -> None:
    _, headers = bootstrap_tenant(client, "Lab Clinical Only")
    patient = client.post("/api/v1/patients", headers=headers, json={}).json()
    sample = client.post(
        "/api/v1/samples",
        headers=headers,
        json={"patient_id": patient["id"]},
    ).json()
    case = client.post(
        "/api/v1/cases",
        headers=headers,
        json={"patient_id": patient["id"], "sample_id": sample["id"]},
    ).json()
    assert case["ml_eligibility"] == "NOT_EVALUATED"
    # Cannot set ML_TRAINING purpose at creation
    bad = client.post(
        "/api/v1/cases",
        headers=headers,
        json={
            "patient_id": patient["id"],
            "sample_id": sample["id"],
            "processing_purposes": ["ML_TRAINING"],
        },
    )
    assert bad.status_code == 422
