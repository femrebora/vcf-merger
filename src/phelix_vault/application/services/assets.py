"""Genomic asset registry, upload sessions, authorized download."""

from __future__ import annotations

import hashlib
import tempfile
from datetime import timedelta
from typing import BinaryIO, cast

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.config import Settings
from phelix_vault.domain.common import (
    AssetStatus,
    AssetType,
    AuditAction,
    DataResidency,
    Permission,
    as_utc,
    new_id,
    utcnow,
)
from phelix_vault.domain.errors import (
    IntegrityError,
    NotFoundError,
    ResidencyViolationError,
    ValidationError,
)
from phelix_vault.infrastructure.database.models import (
    GenomicAssetRow,
    GenomicCaseRow,
    PseudonymousPatientRow,
    SampleRow,
    TenantRow,
    UploadSessionRow,
)
from phelix_vault.infrastructure.kms.base import KeyManagementProvider
from phelix_vault.infrastructure.security import AuthorizationService, Principal
from phelix_vault.infrastructure.storage.base import ObjectStorage, validate_object_key


class AssetService:
    def __init__(
        self,
        session: Session,
        audit: AuditService,
        authz: AuthorizationService,
        storage: ObjectStorage,
        kms: KeyManagementProvider,
        settings: Settings,
    ) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz
        self._storage = storage
        self._kms = kms
        self._settings = settings

    def _bucket(self) -> str:
        if self._storage.provider_name == "s3":
            return self._settings.s3_bucket
        return "vault"

    def create_upload_session(
        self,
        principal: Principal,
        *,
        patient_id: str,
        sample_id: str,
        asset_type: AssetType | str,
        content_type: str,
        genome_build: str | None = None,
        case_id: str | None = None,
        pipeline_name: str | None = None,
        pipeline_version: str | None = None,
        request_id: str | None = None,
    ) -> tuple[GenomicAssetRow, UploadSessionRow]:
        self._authz.require(principal.roles, Permission.ASSET_UPLOAD)
        tenant = self._session.get(TenantRow, principal.tenant_id)
        if not tenant:
            raise NotFoundError("Tenant not found")
        if DataResidency(tenant.required_residency) != self._storage.residency:
            raise ResidencyViolationError(
                "Storage residency does not satisfy tenant data residency policy"
            )

        patient = self._session.get(PseudonymousPatientRow, patient_id)
        sample = self._session.get(SampleRow, sample_id)
        if (
            not patient
            or not sample
            or patient.tenant_id != principal.tenant_id
            or sample.tenant_id != principal.tenant_id
        ):
            raise NotFoundError("Patient or sample not found")
        if case_id:
            case = self._session.get(GenomicCaseRow, case_id)
            if not case or case.tenant_id != principal.tenant_id:
                raise NotFoundError("Case not found")

        atype = asset_type.value if isinstance(asset_type, AssetType) else asset_type
        try:
            AssetType(atype)
        except ValueError as exc:
            raise ValidationError("Invalid asset type") from exc

        asset_id = new_id("PHX-AST")
        # Opaque key — never embed human identifiers.
        object_key = validate_object_key(
            f"tenants/{principal.tenant_id}/assets/{asset_id}/{atype.lower()}"
        )
        data_key = self._kms.create_data_key()
        asset = GenomicAssetRow(
            id=asset_id,
            tenant_id=principal.tenant_id,
            patient_id=patient_id,
            sample_id=sample_id,
            case_id=case_id,
            asset_type=atype,
            storage_provider=self._storage.provider_name,
            storage_bucket=self._bucket(),
            storage_object_key=object_key,
            content_type=content_type,
            genome_build=genome_build,
            pipeline_name=pipeline_name,
            pipeline_version=pipeline_version,
            encryption_key_id=data_key.key_id,
            encrypted_data_key=data_key.encrypted_key,
            encryption_algorithm=data_key.algorithm,
            status=AssetStatus.UPLOAD_PENDING.value,
            residency=self._storage.residency.value,
        )
        self._session.add(asset)
        expires = utcnow() + timedelta(seconds=self._settings.upload_session_ttl_seconds)
        upload_url = self._storage.signed_upload(
            asset.storage_bucket,
            object_key,
            expires_seconds=self._settings.upload_session_ttl_seconds,
        )
        session_row = UploadSessionRow(
            id=new_id("PHX-UPL"),
            tenant_id=principal.tenant_id,
            asset_id=asset.id,
            upload_url=upload_url,
            expires_at=expires,
            created_by=principal.subject,
        )
        self._session.add(session_row)
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.UPLOAD_SESSION_CREATED,
            resource_type="upload_session",
            resource_id=session_row.id,
            request_id=request_id,
            metadata={"asset_id": asset.id, "asset_type": atype},
        )
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.ASSET_CREATED,
            resource_type="asset",
            resource_id=asset.id,
            request_id=request_id,
        )
        return asset, session_row

    def ingest_bytes(
        self,
        principal: Principal,
        session_id: str,
        data: BinaryIO,
        *,
        expected_sha256: str | None = None,
        request_id: str | None = None,
    ) -> GenomicAssetRow:
        """Dev/local finalize path: stream bytes into object storage, verify checksum."""
        self._authz.require(principal.roles, Permission.ASSET_UPLOAD)
        upload = self._session.get(UploadSessionRow, session_id)
        if not upload or upload.tenant_id != principal.tenant_id:
            raise NotFoundError("Upload session not found")
        if upload.completed:
            raise ValidationError("Upload session already completed")
        if as_utc(upload.expires_at) < utcnow():
            raise ValidationError("Upload session expired")

        asset = self._session.get(GenomicAssetRow, upload.asset_id)
        if not asset or asset.tenant_id != principal.tenant_id:
            raise NotFoundError("Asset not found")

        # Hash while buffering to a temp stream — avoid loading multi-GB files as one str.
        # For production S3 direct uploads, bytes arrive via signed URL; this path is for
        # controlled local/dev finalize and modest fixtures.
        hasher = hashlib.sha256()
        total = 0
        with tempfile.SpooledTemporaryFile(max_size=8 * 1024 * 1024) as tmp:
            while True:
                chunk = data.read(1024 * 1024)
                if not chunk:
                    break
                hasher.update(chunk)
                tmp.write(chunk)
                total += len(chunk)
                if total > self._settings.max_upload_bytes:
                    raise ValidationError("Upload exceeds configured maximum size")

            digest = hasher.hexdigest()
            if expected_sha256 and expected_sha256.lower() != digest:
                raise IntegrityError("SHA-256 mismatch")

            tmp.seek(0)
            size = self._storage.put(
                asset.storage_bucket,
                asset.storage_object_key,
                cast(BinaryIO, tmp),
                content_type=asset.content_type,
            )
        verified = self._storage.checksum_sha256(asset.storage_bucket, asset.storage_object_key)
        if verified != digest:
            raise IntegrityError("Stored object checksum mismatch")

        asset.sha256 = digest
        asset.size_bytes = size
        asset.status = AssetStatus.AVAILABLE.value
        asset.ingested_at = utcnow()
        asset.version += 1
        upload.completed = True
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.ASSET_UPLOADED,
            resource_type="asset",
            resource_id=asset.id,
            request_id=request_id,
            metadata={"sha256": digest, "size_bytes": size},
        )
        return asset

    def get(self, principal: Principal, asset_id: str) -> GenomicAssetRow:
        self._authz.require(principal.roles, Permission.ASSET_READ_META)
        asset = self._session.get(GenomicAssetRow, asset_id)
        if not asset or asset.tenant_id != principal.tenant_id:
            # Same message to avoid ID oracle across tenants.
            raise NotFoundError("Asset not found")
        return asset

    def authorize_download(
        self, principal: Principal, asset_id: str, *, request_id: str | None = None
    ) -> dict[str, str | int | None]:
        try:
            self._authz.require(principal.roles, Permission.ASSET_DOWNLOAD)
        except Exception:
            self._audit.record(
                actor_id=principal.subject,
                tenant_id=principal.tenant_id,
                action=AuditAction.ACCESS_DENIED,
                resource_type="asset",
                resource_id=asset_id,
                request_id=request_id,
                result="DENIED",
                metadata={"reason": "missing ASSET_DOWNLOAD"},
            )
            raise

        asset = self._session.get(GenomicAssetRow, asset_id)
        if not asset or asset.tenant_id != principal.tenant_id:
            self._audit.record(
                actor_id=principal.subject,
                tenant_id=principal.tenant_id,
                action=AuditAction.ACCESS_DENIED,
                resource_type="asset",
                resource_id=asset_id,
                request_id=request_id,
                result="DENIED",
                metadata={"reason": "tenant_mismatch_or_missing"},
            )
            raise NotFoundError("Asset not found")
        if asset.status != AssetStatus.AVAILABLE.value:
            raise ValidationError("Asset is not available for download")

        url = self._storage.signed_download(
            asset.storage_bucket,
            asset.storage_object_key,
            expires_seconds=self._settings.signed_url_ttl_seconds,
        )
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.ASSET_DOWNLOADED,
            resource_type="asset",
            resource_id=asset.id,
            request_id=request_id,
            metadata={"expires_seconds": self._settings.signed_url_ttl_seconds},
        )
        return {
            "asset_id": asset.id,
            "download_url": url,
            "expires_seconds": self._settings.signed_url_ttl_seconds,
            "sha256": asset.sha256,
        }

    def mark_pending_deletion(
        self, principal: Principal, asset_id: str, *, request_id: str | None = None
    ) -> GenomicAssetRow:
        self._authz.require(principal.roles, Permission.ASSET_DELETE)
        asset = self.get(principal, asset_id)
        asset.status = AssetStatus.PENDING_DELETION.value
        asset.version += 1
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.ASSET_DELETED,
            resource_type="asset",
            resource_id=asset.id,
            request_id=request_id,
            metadata={"lifecycle": AssetStatus.PENDING_DELETION.value},
        )
        return asset

    def finalize_deletion(
        self, principal: Principal, asset_id: str, *, request_id: str | None = None
    ) -> GenomicAssetRow:
        self._authz.require(principal.roles, Permission.ASSET_DELETE)
        asset = self.get(principal, asset_id)
        if asset.status not in {
            AssetStatus.PENDING_DELETION.value,
            AssetStatus.AVAILABLE.value,
            AssetStatus.UPLOAD_PENDING.value,
        }:
            raise ValidationError("Asset cannot be deleted from current state")
        if self._storage.exists(asset.storage_bucket, asset.storage_object_key):
            self._storage.delete(asset.storage_bucket, asset.storage_object_key)
        asset.status = AssetStatus.DELETED.value
        asset.version += 1
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.ASSET_DELETED,
            resource_type="asset",
            resource_id=asset.id,
            request_id=request_id,
            metadata={
                "lifecycle": AssetStatus.DELETED.value,
                "note": "Logical deletion recorded; physical erasure depends on storage provider",
            },
        )
        return asset

    def list_for_sample(self, principal: Principal, sample_id: str) -> list[GenomicAssetRow]:
        self._authz.require(principal.roles, Permission.ASSET_READ_META)
        return list(
            self._session.scalars(
                select(GenomicAssetRow).where(
                    GenomicAssetRow.tenant_id == principal.tenant_id,
                    GenomicAssetRow.sample_id == sample_id,
                )
            )
        )
