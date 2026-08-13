from __future__ import annotations

from typing import Annotated

from fastapi import APIRouter, Depends, File, Form, UploadFile
from sqlalchemy.orm import Session

from phelix_vault.api.dependencies import build_services, get_db, get_principal, get_request_id
from phelix_vault.api.schemas.common import (
    AssetOut,
    DownloadOut,
    UploadSessionCreate,
    UploadSessionOut,
)
from phelix_vault.infrastructure.security import Principal

router = APIRouter(prefix="/assets", tags=["assets"])


@router.post("/upload-sessions", response_model=UploadSessionOut)
def create_upload_session(
    body: UploadSessionCreate,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> UploadSessionOut:
    services = build_services(session)
    asset, upload = services["assets"].create_upload_session(  # type: ignore[attr-defined]
        principal,
        patient_id=body.patient_id,
        sample_id=body.sample_id,
        asset_type=body.asset_type,
        content_type=body.content_type,
        genome_build=body.genome_build,
        case_id=body.case_id,
        pipeline_name=body.pipeline_name,
        pipeline_version=body.pipeline_version,
        request_id=request_id,
    )
    return UploadSessionOut(
        session_id=upload.id,
        asset_id=asset.id,
        upload_url=upload.upload_url,
        expires_at=upload.expires_at,
        status=asset.status,
    )


@router.post("/upload-sessions/{session_id}/complete", response_model=AssetOut)
async def complete_upload(
    session_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
    file: Annotated[UploadFile, File()],
    expected_sha256: Annotated[str | None, Form()] = None,
) -> AssetOut:
    services = build_services(session)
    # Streaming-friendly: UploadFile.file is a SpooledTemporaryFile.
    asset = services["assets"].ingest_bytes(  # type: ignore[attr-defined]
        principal,
        session_id,
        file.file,
        expected_sha256=expected_sha256,
        request_id=request_id,
    )
    return AssetOut.model_validate(asset)


@router.get("/{asset_id}", response_model=AssetOut)
def get_asset(
    asset_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
) -> AssetOut:
    services = build_services(session)
    asset = services["assets"].get(principal, asset_id)  # type: ignore[attr-defined]
    return AssetOut.model_validate(asset)


@router.post("/{asset_id}/download", response_model=DownloadOut)
def download_asset(
    asset_id: str,
    principal: Annotated[Principal, Depends(get_principal)],
    session: Annotated[Session, Depends(get_db)],
    request_id: Annotated[str, Depends(get_request_id)],
) -> DownloadOut:
    services = build_services(session)
    result = services["assets"].authorize_download(  # type: ignore[attr-defined]
        principal, asset_id, request_id=request_id
    )
    return DownloadOut.model_validate(result)
