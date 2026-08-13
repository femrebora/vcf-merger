from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import IdentityMappingOut, PatientCreate, PatientOut
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/patients", tags=["patients"])


@router.post("", response_model=PatientOut)
def create_patient(
    body: PatientCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> PatientOut:
    services = build_services(session)
    row = services["identity"].create_patient(  # type: ignore[attr-defined]
        principal,
        external_subject_id=body.external_subject_id,
        request_id=request_id,
    )
    return PatientOut.model_validate(row)


@router.get("/{patient_id}", response_model=PatientOut)
def get_patient(
    patient_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> PatientOut:
    services = build_services(session)
    row = services["identity"].get_patient(principal, patient_id)  # type: ignore[attr-defined]
    return PatientOut.model_validate(row)


@router.get("/{patient_id}/identity-mapping", response_model=IdentityMappingOut)
def get_identity_mapping(
    patient_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> IdentityMappingOut:
    services = build_services(session)
    row = services["identity"].get_identity_mapping(  # type: ignore[attr-defined]
        principal, patient_id, request_id=request_id
    )
    return IdentityMappingOut.model_validate(row)
