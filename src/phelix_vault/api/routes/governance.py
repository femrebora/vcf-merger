from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal
from phelix_vault.api.schemas.common import AuthorizationCreate, AuthorizationOut
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/governance", tags=["governance"])


@router.post("/authorizations", response_model=AuthorizationOut)
def create_authorization(
    body: AuthorizationCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> AuthorizationOut:
    services = build_services(session)
    row = services["governance"].create_authorization(  # type: ignore[attr-defined]
        principal,
        case_id=body.case_id,
        processing_purpose=body.processing_purpose,
        legal_basis_code=body.legal_basis_code,
        authorization_source=body.authorization_source,
        international_transfer_allowed=body.international_transfer_allowed,
        ml_training_allowed=body.ml_training_allowed,
        research_allowed=body.research_allowed,
        reanalysis_allowed=body.reanalysis_allowed,
        notes=body.notes,
    )
    return AuthorizationOut.model_validate(row)
