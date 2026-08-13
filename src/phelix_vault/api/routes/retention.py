from __future__ import annotations

from typing import Annotated, Any

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import RetentionPolicyCreate, RetentionPolicyOut
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/retention", tags=["retention"])


@router.post("/policies", response_model=RetentionPolicyOut)
def create_policy(
    body: RetentionPolicyCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> RetentionPolicyOut:
    services = build_services(session)
    row = services["retention"].create_policy(  # type: ignore[attr-defined]
        principal,
        name=body.name,
        asset_type=body.asset_type,
        retention_period_days=body.retention_period_days,
        action=body.action,
        legal_hold_supported=body.legal_hold_supported,
        request_id=request_id,
    )
    return RetentionPolicyOut.model_validate(row)


@router.get("/evaluate")
def evaluate(
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> list[dict[str, Any]]:
    services = build_services(session)
    due = services["retention"].evaluate(principal)  # type: ignore[attr-defined]
    return list(due)
