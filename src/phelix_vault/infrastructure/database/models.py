"""Relational models — domain boundaries kept as separate tables."""

from __future__ import annotations

from datetime import datetime
from typing import Any

from sqlalchemy import (
    JSON,
    Boolean,
    DateTime,
    ForeignKey,
    Integer,
    LargeBinary,
    String,
    Text,
    UniqueConstraint,
)
from sqlalchemy.orm import Mapped, mapped_column

from phelix_vault.domain.common import utcnow
from phelix_vault.infrastructure.database.base import Base


class TenantRow(Base):
    __tablename__ = "tenants"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    required_residency: Mapped[str] = mapped_column(String(16), nullable=False, default="TR")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )
    version: Mapped[int] = mapped_column(Integer, default=1, nullable=False)


class PseudonymousPatientRow(Base):
    __tablename__ = "pseudonymous_patients"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )
    version: Mapped[int] = mapped_column(Integer, default=1, nullable=False)


class IdentityMappingRow(Base):
    """Isolated identity domain — not joined casually to genomic asset queries."""

    __tablename__ = "identity_mappings"
    __table_args__ = (UniqueConstraint("tenant_id", "external_subject_id", name="uq_identity_ext"),)

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    patient_id: Mapped[str] = mapped_column(
        ForeignKey("pseudonymous_patients.id"), unique=True, nullable=False
    )
    # Encrypted / opaque external reference — never use as object storage filename.
    external_subject_id: Mapped[str] = mapped_column(String(512), nullable=False)
    display_label_encrypted: Mapped[bytes | None] = mapped_column(LargeBinary, nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )


class SampleRow(Base):
    __tablename__ = "samples"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    patient_id: Mapped[str] = mapped_column(
        ForeignKey("pseudonymous_patients.id"), index=True, nullable=False
    )
    sample_type: Mapped[str] = mapped_column(String(64), default="WES")
    collection_date: Mapped[str | None] = mapped_column(String(32), nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )
    version: Mapped[int] = mapped_column(Integer, default=1, nullable=False)


class GenomicCaseRow(Base):
    __tablename__ = "genomic_cases"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    patient_id: Mapped[str] = mapped_column(
        ForeignKey("pseudonymous_patients.id"), index=True, nullable=False
    )
    sample_id: Mapped[str] = mapped_column(ForeignKey("samples.id"), index=True, nullable=False)
    curation_status: Mapped[str] = mapped_column(String(64), default="RAW")
    ml_eligibility: Mapped[str] = mapped_column(String(64), default="NOT_EVALUATED")
    data_zone: Mapped[str] = mapped_column(String(64), default="RAW_CLINICAL")
    processing_purposes: Mapped[list[Any]] = mapped_column(JSON, default=list)
    # Data quality dimensions (explicit, not inferred silently)
    technical_completeness: Mapped[bool] = mapped_column(Boolean, default=False)
    phenotype_completeness: Mapped[bool] = mapped_column(Boolean, default=False)
    report_available: Mapped[bool] = mapped_column(Boolean, default=False)
    diagnosis_available: Mapped[bool] = mapped_column(Boolean, default=False)
    variant_confirmed: Mapped[bool] = mapped_column(Boolean, default=False)
    expert_review_count: Mapped[int] = mapped_column(Integer, default=0)
    pipeline_provenance_known: Mapped[bool] = mapped_column(Boolean, default=False)
    # Technical vs clinical truth are separate metadata slots
    technical_truth_notes: Mapped[str | None] = mapped_column(Text, nullable=True)
    clinical_truth_notes: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )
    version: Mapped[int] = mapped_column(Integer, default=1, nullable=False)


class CaseStatusTransitionRow(Base):
    __tablename__ = "case_status_transitions"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    case_id: Mapped[str] = mapped_column(ForeignKey("genomic_cases.id"), index=True, nullable=False)
    field_name: Mapped[str] = mapped_column(String(64), nullable=False)
    previous_state: Mapped[str] = mapped_column(String(64), nullable=False)
    new_state: Mapped[str] = mapped_column(String(64), nullable=False)
    reason: Mapped[str | None] = mapped_column(Text, nullable=True)
    actor_id: Mapped[str] = mapped_column(String(128), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class RetentionPolicyRow(Base):
    __tablename__ = "retention_policies"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    asset_type: Mapped[str] = mapped_column(String(64), nullable=False)
    retention_period_days: Mapped[int] = mapped_column(Integer, nullable=False)
    action: Mapped[str] = mapped_column(String(32), nullable=False)
    legal_hold_supported: Mapped[bool] = mapped_column(Boolean, default=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class GenomicAssetRow(Base):
    __tablename__ = "genomic_assets"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    patient_id: Mapped[str] = mapped_column(
        ForeignKey("pseudonymous_patients.id"), index=True, nullable=False
    )
    sample_id: Mapped[str] = mapped_column(ForeignKey("samples.id"), index=True, nullable=False)
    case_id: Mapped[str | None] = mapped_column(
        ForeignKey("genomic_cases.id"), index=True, nullable=True
    )
    asset_type: Mapped[str] = mapped_column(String(64), nullable=False)
    storage_provider: Mapped[str] = mapped_column(String(64), nullable=False)
    storage_bucket: Mapped[str] = mapped_column(String(255), nullable=False)
    storage_object_key: Mapped[str] = mapped_column(String(1024), nullable=False)
    size_bytes: Mapped[int | None] = mapped_column(Integer, nullable=True)
    sha256: Mapped[str | None] = mapped_column(String(64), nullable=True)
    content_type: Mapped[str | None] = mapped_column(String(128), nullable=True)
    genome_build: Mapped[str | None] = mapped_column(String(64), nullable=True)
    pipeline_name: Mapped[str | None] = mapped_column(String(128), nullable=True)
    pipeline_version: Mapped[str | None] = mapped_column(String(64), nullable=True)
    encryption_key_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    encrypted_data_key: Mapped[bytes | None] = mapped_column(LargeBinary, nullable=True)
    encryption_algorithm: Mapped[str | None] = mapped_column(String(64), nullable=True)
    retention_policy_id: Mapped[str | None] = mapped_column(
        ForeignKey("retention_policies.id"), nullable=True
    )
    status: Mapped[str] = mapped_column(String(64), default="UPLOAD_PENDING")
    data_zone: Mapped[str] = mapped_column(String(64), default="RAW_CLINICAL")
    residency: Mapped[str] = mapped_column(String(16), default="TR")
    relation_kind: Mapped[str] = mapped_column(String(32), default="ORIGINAL")
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    ingested_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )
    version: Mapped[int] = mapped_column(Integer, default=1, nullable=False)


class AssetProvenanceRow(Base):
    __tablename__ = "asset_provenance"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    parent_asset_id: Mapped[str] = mapped_column(
        ForeignKey("genomic_assets.id"), index=True, nullable=False
    )
    child_asset_id: Mapped[str] = mapped_column(
        ForeignKey("genomic_assets.id"), index=True, nullable=False
    )
    relation_kind: Mapped[str] = mapped_column(String(32), nullable=False)
    pipeline_run_id: Mapped[str | None] = mapped_column(String(128), nullable=True)
    workflow_engine: Mapped[str | None] = mapped_column(String(128), nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class UploadSessionRow(Base):
    __tablename__ = "upload_sessions"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    asset_id: Mapped[str] = mapped_column(ForeignKey("genomic_assets.id"), nullable=False)
    upload_url: Mapped[str] = mapped_column(Text, nullable=False)
    expires_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)
    completed: Mapped[bool] = mapped_column(Boolean, default=False)
    created_by: Mapped[str] = mapped_column(String(128), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class ProcessingAuthorizationRow(Base):
    __tablename__ = "processing_authorizations"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    case_id: Mapped[str] = mapped_column(ForeignKey("genomic_cases.id"), index=True, nullable=False)
    processing_purpose: Mapped[str] = mapped_column(String(64), nullable=False)
    legal_basis_code: Mapped[str] = mapped_column(String(128), nullable=False)
    authorization_source: Mapped[str] = mapped_column(String(255), nullable=False)
    effective_from: Mapped[datetime] = mapped_column(DateTime(timezone=True), nullable=False)
    effective_until: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    restrictions: Mapped[str | None] = mapped_column(Text, nullable=True)
    data_controller: Mapped[str | None] = mapped_column(String(255), nullable=True)
    data_processor: Mapped[str | None] = mapped_column(String(255), nullable=True)
    international_transfer_allowed: Mapped[bool] = mapped_column(Boolean, default=False)
    ml_training_allowed: Mapped[bool] = mapped_column(Boolean, default=False)
    research_allowed: Mapped[bool] = mapped_column(Boolean, default=False)
    reanalysis_allowed: Mapped[bool] = mapped_column(Boolean, default=True)
    approved_by: Mapped[str] = mapped_column(String(128), nullable=False)
    approved_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    notes: Mapped[str | None] = mapped_column(Text, nullable=True)


class DatasetRow(Base):
    __tablename__ = "datasets"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    dataset_role: Mapped[str] = mapped_column(String(64), nullable=False)
    source: Mapped[str] = mapped_column(String(64), default="INTERNAL")
    created_by: Mapped[str] = mapped_column(String(128), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)
    updated_at: Mapped[datetime] = mapped_column(
        DateTime(timezone=True), default=utcnow, onupdate=utcnow
    )


class DatasetVersionRow(Base):
    __tablename__ = "dataset_versions"
    __table_args__ = (UniqueConstraint("dataset_id", "version", name="uq_dataset_version"),)

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    dataset_id: Mapped[str] = mapped_column(ForeignKey("datasets.id"), index=True, nullable=False)
    version: Mapped[str] = mapped_column(String(32), nullable=False)
    purpose: Mapped[str] = mapped_column(String(64), nullable=False)
    dataset_role: Mapped[str] = mapped_column(String(64), nullable=False)
    approval_state: Mapped[str] = mapped_column(String(64), default="DRAFT")
    eligibility_criteria: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    excluded_case_ids: Mapped[list[Any]] = mapped_column(JSON, default=list)
    created_by: Mapped[str] = mapped_column(String(128), nullable=False)
    approved_by: Mapped[str | None] = mapped_column(String(128), nullable=True)
    approved_at: Mapped[datetime | None] = mapped_column(DateTime(timezone=True), nullable=True)
    immutable: Mapped[bool] = mapped_column(Boolean, default=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class DatasetCaseRow(Base):
    __tablename__ = "dataset_cases"
    __table_args__ = (UniqueConstraint("dataset_version_id", "case_id", name="uq_dataset_case"),)

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    dataset_version_id: Mapped[str] = mapped_column(
        ForeignKey("dataset_versions.id"), index=True, nullable=False
    )
    case_id: Mapped[str] = mapped_column(ForeignKey("genomic_cases.id"), index=True, nullable=False)
    sample_id: Mapped[str] = mapped_column(ForeignKey("samples.id"), nullable=False)
    patient_id: Mapped[str] = mapped_column(ForeignKey("pseudonymous_patients.id"), nullable=False)
    vcf_asset_id: Mapped[str | None] = mapped_column(ForeignKey("genomic_assets.id"), nullable=True)
    phenotype_ref: Mapped[str | None] = mapped_column(String(255), nullable=True)
    truth_label_ref: Mapped[str | None] = mapped_column(String(255), nullable=True)


class ExternalCohortRow(Base):
    __tablename__ = "external_cohorts"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    name: Mapped[str] = mapped_column(String(255), nullable=False)
    description: Mapped[str | None] = mapped_column(Text, nullable=True)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class ExternalCohortVersionRow(Base):
    __tablename__ = "external_cohort_versions"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    cohort_id: Mapped[str] = mapped_column(ForeignKey("external_cohorts.id"), index=True)
    version: Mapped[str] = mapped_column(String(32), nullable=False)
    dataset_version_id: Mapped[str | None] = mapped_column(
        ForeignKey("dataset_versions.id"), nullable=True
    )
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class ModelDatasetLinkRow(Base):
    """Generic reference: DatasetVersion used by ModelVersion (no ML deps)."""

    __tablename__ = "model_dataset_links"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    tenant_id: Mapped[str] = mapped_column(ForeignKey("tenants.id"), index=True, nullable=False)
    dataset_version_id: Mapped[str] = mapped_column(
        ForeignKey("dataset_versions.id"), nullable=False
    )
    model_version_ref: Mapped[str] = mapped_column(String(255), nullable=False)
    created_at: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow)


class AuditEventRow(Base):
    __tablename__ = "audit_events"

    id: Mapped[str] = mapped_column(String(64), primary_key=True)
    timestamp: Mapped[datetime] = mapped_column(DateTime(timezone=True), default=utcnow, index=True)
    actor_id: Mapped[str] = mapped_column(String(128), nullable=False)
    tenant_id: Mapped[str] = mapped_column(String(64), index=True, nullable=False)
    action: Mapped[str] = mapped_column(String(64), nullable=False)
    resource_type: Mapped[str] = mapped_column(String(64), nullable=False)
    resource_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    request_id: Mapped[str | None] = mapped_column(String(64), nullable=True)
    source_ip: Mapped[str | None] = mapped_column(String(64), nullable=True)
    result: Mapped[str] = mapped_column(String(32), nullable=False)
    metadata_json: Mapped[dict[str, Any]] = mapped_column(JSON, default=dict)
    prev_hash: Mapped[str | None] = mapped_column(String(64), nullable=True)
    event_hash: Mapped[str] = mapped_column(String(64), nullable=False)
