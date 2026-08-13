from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import SampleCreate, SampleOut
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/samples", tags=["samples"])


@router.post("", response_model=SampleOut)
def create_sample(
    body: SampleCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> SampleOut:
    services = build_services(session)
    row = services["samples"].create(  # type: ignore[attr-defined]
        principal,
        patient_id=body.patient_id,
        sample_type=body.sample_type,
        collection_date=body.collection_date,
        request_id=request_id,
    )
    return SampleOut.model_validate(row)


@router.get("/{sample_id}", response_model=SampleOut)
def get_sample(
    sample_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> SampleOut:
    services = build_services(session)
    row = services["samples"].get(principal, sample_id)  # type: ignore[attr-defined]
    return SampleOut.model_validate(row)
