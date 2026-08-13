from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal
from phelix_vault.api.schemas.common import AuditEventOut
from phelix_vault.domain.common import Permission
from phelix_vault.infrastructure.security import AuthorizationService, Principal

router = APIRouter(prefix="/audit", tags=["audit"])


@router.get("/events", response_model=list[AuditEventOut])
def list_events(
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> list[AuditEventOut]:
    AuthorizationService().require(principal.roles, Permission.AUDIT_READ)
    services = build_services(session)
    rows = services["audit"].list_for_tenant(principal.tenant_id)  # type: ignore[attr-defined]
    return [AuditEventOut.model_validate(r) for r in rows]
