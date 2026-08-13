"""Development-only token minting. Disabled conceptually in production deployments."""

from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends, HTTPException

from phelix_vault.api.dependencies import get_container
from phelix_vault.api.schemas.common import DevTokenRequest, DevTokenResponse
from phelix_vault.config import Settings, get_settings

router = APIRouter(prefix="/auth", tags=["auth"])


@router.post("/dev-token", response_model=DevTokenResponse)
def mint_dev_token(
    body: DevTokenRequest,
    settings: Annotated[Settings, Depends(get_settings)],
) -> DevTokenResponse:
    if settings.environment == "production":
        raise HTTPException(status_code=404, detail="Not found")
    container = get_container()
    token = container.token_service.create_access_token(
        subject=body.subject,
        tenant_id=body.tenant_id,
        roles=body.roles,
        email=body.email,
    )
    return DevTokenResponse(access_token=token)
