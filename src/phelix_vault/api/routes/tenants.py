from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import TenantCreate, TenantOut
from phelix_vault.domain.common import Permission
from phelix_vault.infrastructure.security import AuthorizationService, Principal

router = APIRouter(prefix="/tenants", tags=["tenants"])


@router.post("", response_model=TenantOut)
def create_tenant(
    body: TenantCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> TenantOut:
    # Bootstrap: SYSTEM_ADMIN or TENANT_ADMIN may create; first tenant may use any authed admin.
    AuthorizationService().require(principal.roles, Permission.TENANT_MANAGE)
    services = build_services(session)
    row = services["tenants"].create(  # type: ignore[attr-defined]
        name=body.name,
        required_residency=body.required_residency,
        actor_id=principal.subject,
        request_id=request_id,
    )
    return TenantOut.model_validate(row)


@router.get("/{tenant_id}", response_model=TenantOut)
def get_tenant(
    tenant_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> TenantOut:
    if principal.tenant_id != tenant_id:
        from phelix_vault.domain.errors import AuthorizationError

        raise AuthorizationError("Cross-tenant access denied")
    services = build_services(session)
    row = services["tenants"].get(tenant_id)  # type: ignore[attr-defined]
    return TenantOut.model_validate(row)
