from __future__ import annotations

from typing import Annotated, Any

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import (
    DatasetCreate,
    DatasetOut,
    DatasetVersionCreate,
    DatasetVersionOut,
)
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/datasets", tags=["datasets"])


@router.post("", response_model=DatasetOut)
def create_dataset(
    body: DatasetCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> DatasetOut:
    services = build_services(session)
    row = services["datasets"].create_dataset(  # type: ignore[attr-defined]
        principal,
        name=body.name,
        dataset_role=body.dataset_role,
        source=body.source,
        request_id=request_id,
    )
    return DatasetOut.model_validate(row)


@router.post("/{dataset_id}/versions", response_model=DatasetVersionOut)
def create_version(
    dataset_id: str,
    body: DatasetVersionCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> DatasetVersionOut:
    services = build_services(session)
    row = services["datasets"].create_version(  # type: ignore[attr-defined]
        principal,
        dataset_id,
        version=body.version,
        purpose=body.purpose,
        case_ids=body.case_ids,
        excluded_case_ids=body.excluded_case_ids,
        eligibility_criteria=body.eligibility_criteria,
        request_id=request_id,
    )
    return DatasetVersionOut.model_validate(row)


@router.post("/versions/{dataset_version_id}/approve", response_model=DatasetVersionOut)
def approve_version(
    dataset_version_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> DatasetVersionOut:
    services = build_services(session)
    row = services["datasets"].approve_version(  # type: ignore[attr-defined]
        principal, dataset_version_id, request_id=request_id
    )
    return DatasetVersionOut.model_validate(row)


@router.get("/versions/{dataset_version_id}/manifest")
def export_manifest(
    dataset_version_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> dict[str, Any]:
    services = build_services(session)
    manifest = services["datasets"].export_manifest(  # type: ignore[attr-defined]
        principal, dataset_version_id, request_id=request_id
    )
    return dict(manifest)


@router.get("/versions/{version_a}/overlap/{version_b}")
def overlap(
    version_a: str,
    version_b: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> dict[str, Any]:
    services = build_services(session)
    report = services["datasets"].patient_overlap_report(  # type: ignore[attr-defined]
        principal, version_a, version_b
    )
    return dict(report)
