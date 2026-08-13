"""API DTOs — never expose ORM models directly."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from pydantic import BaseModel, ConfigDict, Field


class ORMModel(BaseModel):
    model_config = ConfigDict(from_attributes=True)


class TenantCreate(BaseModel):
    name: str = Field(min_length=1, max_length=255)
    required_residency: str = "TR"


class TenantOut(ORMModel):
    id: str
    name: str
    required_residency: str
    created_at: datetime
    version: int


class PatientCreate(BaseModel):
    external_subject_id: str | None = None


class PatientOut(ORMModel):
    id: str
    tenant_id: str
    created_at: datetime


class IdentityMappingOut(ORMModel):
    id: str
    patient_id: str
    external_subject_id: str
    created_at: datetime


class SampleCreate(BaseModel):
    patient_id: str
    sample_type: str = "WES"
    collection_date: str | None = None


class SampleOut(ORMModel):
    id: str
    tenant_id: str
    patient_id: str
    sample_type: str
    collection_date: str | None
    created_at: datetime


class CaseCreate(BaseModel):
    patient_id: str
    sample_id: str
    processing_purposes: list[str] = Field(default_factory=lambda: ["CLINICAL_DIAGNOSIS"])


class CaseOut(ORMModel):
    id: str
    tenant_id: str
    patient_id: str
    sample_id: str
    curation_status: str
    ml_eligibility: str
    data_zone: str
    processing_purposes: list[Any]
    version: int


class CurationUpdate(BaseModel):
    status: str
    reason: str | None = None


class MlEligibilityUpdate(BaseModel):
    status: str
    reason: str | None = None


class UploadSessionCreate(BaseModel):
    patient_id: str
    sample_id: str
    asset_type: str
    content_type: str = "application/octet-stream"
    genome_build: str | None = None
    case_id: str | None = None
    pipeline_name: str | None = None
    pipeline_version: str | None = None


class UploadSessionOut(BaseModel):
    session_id: str
    asset_id: str
    upload_url: str
    expires_at: datetime
    status: str


class AssetOut(ORMModel):
    id: str
    tenant_id: str
    patient_id: str
    sample_id: str
    case_id: str | None
    asset_type: str
    storage_provider: str
    storage_bucket: str
    storage_object_key: str
    size_bytes: int | None
    sha256: str | None
    content_type: str | None
    genome_build: str | None
    status: str
    residency: str
    created_at: datetime
    ingested_at: datetime | None


class DownloadOut(BaseModel):
    asset_id: str
    download_url: str
    expires_seconds: int
    sha256: str | None


class DatasetCreate(BaseModel):
    name: str
    dataset_role: str
    source: str = "INTERNAL"


class DatasetOut(ORMModel):
    id: str
    tenant_id: str
    name: str
    dataset_role: str
    source: str
    created_by: str
    created_at: datetime


class DatasetVersionCreate(BaseModel):
    version: str
    purpose: str
    case_ids: list[str]
    excluded_case_ids: list[str] = Field(default_factory=list)
    eligibility_criteria: dict[str, Any] = Field(default_factory=dict)


class DatasetVersionOut(ORMModel):
    id: str
    dataset_id: str
    version: str
    purpose: str
    dataset_role: str
    approval_state: str
    created_by: str
    approved_by: str | None
    immutable: bool
    created_at: datetime


class AuthorizationCreate(BaseModel):
    case_id: str
    processing_purpose: str
    legal_basis_code: str
    authorization_source: str
    international_transfer_allowed: bool = False
    ml_training_allowed: bool = False
    research_allowed: bool = False
    reanalysis_allowed: bool = True
    notes: str | None = None


class AuthorizationOut(ORMModel):
    id: str
    case_id: str
    processing_purpose: str
    legal_basis_code: str
    ml_training_allowed: bool
    approved_by: str
    approved_at: datetime


class RetentionPolicyCreate(BaseModel):
    name: str
    asset_type: str
    retention_period_days: int
    action: str
    legal_hold_supported: bool = True


class RetentionPolicyOut(ORMModel):
    id: str
    name: str
    asset_type: str
    retention_period_days: int
    action: str


class AuditEventOut(ORMModel):
    id: str
    timestamp: datetime
    actor_id: str
    tenant_id: str
    action: str
    resource_type: str
    resource_id: str | None
    result: str
    event_hash: str


class DevTokenRequest(BaseModel):
    subject: str = "dev-user"
    tenant_id: str
    roles: list[str]
    email: str | None = None


class DevTokenResponse(BaseModel):
    access_token: str
    token_type: str = "bearer"


class ErrorOut(BaseModel):
    error: str
    message: str
    request_id: str | None = None
