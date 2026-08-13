from __future__ import annotations

from typing import Annotated, Any

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import CaseCreate, CaseOut, CurationUpdate, MlEligibilityUpdate
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/cases", tags=["cases"])


@router.post("", response_model=CaseOut)
def create_case(
    body: CaseCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> CaseOut:
    services = build_services(session)
    row = services["cases"].create(  # type: ignore[attr-defined]
        principal,
        patient_id=body.patient_id,
        sample_id=body.sample_id,
        purposes=body.processing_purposes,
        request_id=request_id,
    )
    return CaseOut.model_validate(row)


@router.get("/{case_id}", response_model=CaseOut)
def get_case(
    case_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> CaseOut:
    services = build_services(session)
    row = services["cases"].get(principal, case_id)  # type: ignore[attr-defined]
    return CaseOut.model_validate(row)


@router.post("/{case_id}/curation-status", response_model=CaseOut)
def set_curation(
    case_id: str,
    body: CurationUpdate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> CaseOut:
    services = build_services(session)
    row = services["cases"].set_curation_status(  # type: ignore[attr-defined]
        principal, case_id, status=body.status, reason=body.reason, request_id=request_id
    )
    return CaseOut.model_validate(row)


@router.post("/{case_id}/ml-eligibility", response_model=CaseOut)
def set_ml_eligibility(
    case_id: str,
    body: MlEligibilityUpdate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> CaseOut:
    services = build_services(session)
    row = services["cases"].set_ml_eligibility(  # type: ignore[attr-defined]
        principal, case_id, status=body.status, reason=body.reason, request_id=request_id
    )
    return CaseOut.model_validate(row)


@router.get("/{case_id}/eligibility-report")
def eligibility_report(
    case_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> dict[str, Any]:
    services = build_services(session)
    report = services["cases"].eligibility_report(principal, case_id)  # type: ignore[attr-defined]
    return dict(report)
