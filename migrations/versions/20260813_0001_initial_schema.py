"""Initial pHelix Vault schema.

Revision ID: 20260813_0001
Revises:
Create Date: 2026-08-13

"""

from __future__ import annotations

from typing import Sequence, Union

from alembic import op
import sqlalchemy as sa

revision: str = "20260813_0001"
down_revision: Union[str, None] = None
branch_labels: Union[str, Sequence[str], None] = None
depends_on: Union[str, Sequence[str], None] = None


def upgrade() -> None:
    op.create_table(
        "tenants",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("required_residency", sa.String(length=16), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.Column("version", sa.Integer(), nullable=False),
    )
    op.create_table(
        "pseudonymous_patients",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.Column("version", sa.Integer(), nullable=False),
    )
    op.create_index("ix_pseudonymous_patients_tenant_id", "pseudonymous_patients", ["tenant_id"])
    op.create_table(
        "identity_mappings",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "patient_id",
            sa.String(length=64),
            sa.ForeignKey("pseudonymous_patients.id"),
            nullable=False,
            unique=True,
        ),
        sa.Column("external_subject_id", sa.String(length=512), nullable=False),
        sa.Column("display_label_encrypted", sa.LargeBinary(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.UniqueConstraint("tenant_id", "external_subject_id", name="uq_identity_ext"),
    )
    op.create_table(
        "samples",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "patient_id",
            sa.String(length=64),
            sa.ForeignKey("pseudonymous_patients.id"),
            nullable=False,
        ),
        sa.Column("sample_type", sa.String(length=64)),
        sa.Column("collection_date", sa.String(length=32), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.Column("version", sa.Integer(), nullable=False),
    )
    op.create_table(
        "retention_policies",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("asset_type", sa.String(length=64), nullable=False),
        sa.Column("retention_period_days", sa.Integer(), nullable=False),
        sa.Column("action", sa.String(length=32), nullable=False),
        sa.Column("legal_hold_supported", sa.Boolean()),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "genomic_cases",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "patient_id",
            sa.String(length=64),
            sa.ForeignKey("pseudonymous_patients.id"),
            nullable=False,
        ),
        sa.Column("sample_id", sa.String(length=64), sa.ForeignKey("samples.id"), nullable=False),
        sa.Column("curation_status", sa.String(length=64)),
        sa.Column("ml_eligibility", sa.String(length=64)),
        sa.Column("data_zone", sa.String(length=64)),
        sa.Column("processing_purposes", sa.JSON()),
        sa.Column("technical_completeness", sa.Boolean()),
        sa.Column("phenotype_completeness", sa.Boolean()),
        sa.Column("report_available", sa.Boolean()),
        sa.Column("diagnosis_available", sa.Boolean()),
        sa.Column("variant_confirmed", sa.Boolean()),
        sa.Column("expert_review_count", sa.Integer()),
        sa.Column("pipeline_provenance_known", sa.Boolean()),
        sa.Column("technical_truth_notes", sa.Text(), nullable=True),
        sa.Column("clinical_truth_notes", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.Column("version", sa.Integer(), nullable=False),
    )
    op.create_table(
        "case_status_transitions",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("case_id", sa.String(length=64), sa.ForeignKey("genomic_cases.id"), nullable=False),
        sa.Column("field_name", sa.String(length=64), nullable=False),
        sa.Column("previous_state", sa.String(length=64), nullable=False),
        sa.Column("new_state", sa.String(length=64), nullable=False),
        sa.Column("reason", sa.Text(), nullable=True),
        sa.Column("actor_id", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "genomic_assets",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "patient_id",
            sa.String(length=64),
            sa.ForeignKey("pseudonymous_patients.id"),
            nullable=False,
        ),
        sa.Column("sample_id", sa.String(length=64), sa.ForeignKey("samples.id"), nullable=False),
        sa.Column("case_id", sa.String(length=64), sa.ForeignKey("genomic_cases.id"), nullable=True),
        sa.Column("asset_type", sa.String(length=64), nullable=False),
        sa.Column("storage_provider", sa.String(length=64), nullable=False),
        sa.Column("storage_bucket", sa.String(length=255), nullable=False),
        sa.Column("storage_object_key", sa.String(length=1024), nullable=False),
        sa.Column("size_bytes", sa.Integer(), nullable=True),
        sa.Column("sha256", sa.String(length=64), nullable=True),
        sa.Column("content_type", sa.String(length=128), nullable=True),
        sa.Column("genome_build", sa.String(length=64), nullable=True),
        sa.Column("pipeline_name", sa.String(length=128), nullable=True),
        sa.Column("pipeline_version", sa.String(length=64), nullable=True),
        sa.Column("encryption_key_id", sa.String(length=128), nullable=True),
        sa.Column("encrypted_data_key", sa.LargeBinary(), nullable=True),
        sa.Column("encryption_algorithm", sa.String(length=64), nullable=True),
        sa.Column(
            "retention_policy_id",
            sa.String(length=64),
            sa.ForeignKey("retention_policies.id"),
            nullable=True,
        ),
        sa.Column("status", sa.String(length=64)),
        sa.Column("data_zone", sa.String(length=64)),
        sa.Column("residency", sa.String(length=16)),
        sa.Column("relation_kind", sa.String(length=32)),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("ingested_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
        sa.Column("version", sa.Integer(), nullable=False),
    )
    op.create_table(
        "asset_provenance",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "parent_asset_id", sa.String(length=64), sa.ForeignKey("genomic_assets.id"), nullable=False
        ),
        sa.Column(
            "child_asset_id", sa.String(length=64), sa.ForeignKey("genomic_assets.id"), nullable=False
        ),
        sa.Column("relation_kind", sa.String(length=32), nullable=False),
        sa.Column("pipeline_run_id", sa.String(length=128), nullable=True),
        sa.Column("workflow_engine", sa.String(length=128), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "upload_sessions",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("asset_id", sa.String(length=64), sa.ForeignKey("genomic_assets.id"), nullable=False),
        sa.Column("upload_url", sa.Text(), nullable=False),
        sa.Column("expires_at", sa.DateTime(timezone=True), nullable=False),
        sa.Column("completed", sa.Boolean()),
        sa.Column("created_by", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "processing_authorizations",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("case_id", sa.String(length=64), sa.ForeignKey("genomic_cases.id"), nullable=False),
        sa.Column("processing_purpose", sa.String(length=64), nullable=False),
        sa.Column("legal_basis_code", sa.String(length=128), nullable=False),
        sa.Column("authorization_source", sa.String(length=255), nullable=False),
        sa.Column("effective_from", sa.DateTime(timezone=True), nullable=False),
        sa.Column("effective_until", sa.DateTime(timezone=True), nullable=True),
        sa.Column("restrictions", sa.Text(), nullable=True),
        sa.Column("data_controller", sa.String(length=255), nullable=True),
        sa.Column("data_processor", sa.String(length=255), nullable=True),
        sa.Column("international_transfer_allowed", sa.Boolean()),
        sa.Column("ml_training_allowed", sa.Boolean()),
        sa.Column("research_allowed", sa.Boolean()),
        sa.Column("reanalysis_allowed", sa.Boolean()),
        sa.Column("approved_by", sa.String(length=128), nullable=False),
        sa.Column("approved_at", sa.DateTime(timezone=True)),
        sa.Column("notes", sa.Text(), nullable=True),
    )
    op.create_table(
        "datasets",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("dataset_role", sa.String(length=64), nullable=False),
        sa.Column("source", sa.String(length=64)),
        sa.Column("created_by", sa.String(length=128), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.Column("updated_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "dataset_versions",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("dataset_id", sa.String(length=64), sa.ForeignKey("datasets.id"), nullable=False),
        sa.Column("version", sa.String(length=32), nullable=False),
        sa.Column("purpose", sa.String(length=64), nullable=False),
        sa.Column("dataset_role", sa.String(length=64), nullable=False),
        sa.Column("approval_state", sa.String(length=64)),
        sa.Column("eligibility_criteria", sa.JSON()),
        sa.Column("excluded_case_ids", sa.JSON()),
        sa.Column("created_by", sa.String(length=128), nullable=False),
        sa.Column("approved_by", sa.String(length=128), nullable=True),
        sa.Column("approved_at", sa.DateTime(timezone=True), nullable=True),
        sa.Column("immutable", sa.Boolean()),
        sa.Column("created_at", sa.DateTime(timezone=True)),
        sa.UniqueConstraint("dataset_id", "version", name="uq_dataset_version"),
    )
    op.create_table(
        "dataset_cases",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "dataset_version_id",
            sa.String(length=64),
            sa.ForeignKey("dataset_versions.id"),
            nullable=False,
        ),
        sa.Column("case_id", sa.String(length=64), sa.ForeignKey("genomic_cases.id"), nullable=False),
        sa.Column("sample_id", sa.String(length=64), sa.ForeignKey("samples.id"), nullable=False),
        sa.Column(
            "patient_id",
            sa.String(length=64),
            sa.ForeignKey("pseudonymous_patients.id"),
            nullable=False,
        ),
        sa.Column(
            "vcf_asset_id", sa.String(length=64), sa.ForeignKey("genomic_assets.id"), nullable=True
        ),
        sa.Column("phenotype_ref", sa.String(length=255), nullable=True),
        sa.Column("truth_label_ref", sa.String(length=255), nullable=True),
        sa.UniqueConstraint("dataset_version_id", "case_id", name="uq_dataset_case"),
    )
    op.create_table(
        "external_cohorts",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("name", sa.String(length=255), nullable=False),
        sa.Column("description", sa.Text(), nullable=True),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "external_cohort_versions",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column("cohort_id", sa.String(length=64), sa.ForeignKey("external_cohorts.id")),
        sa.Column("version", sa.String(length=32), nullable=False),
        sa.Column(
            "dataset_version_id",
            sa.String(length=64),
            sa.ForeignKey("dataset_versions.id"),
            nullable=True,
        ),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "model_dataset_links",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("tenant_id", sa.String(length=64), sa.ForeignKey("tenants.id"), nullable=False),
        sa.Column(
            "dataset_version_id",
            sa.String(length=64),
            sa.ForeignKey("dataset_versions.id"),
            nullable=False,
        ),
        sa.Column("model_version_ref", sa.String(length=255), nullable=False),
        sa.Column("created_at", sa.DateTime(timezone=True)),
    )
    op.create_table(
        "audit_events",
        sa.Column("id", sa.String(length=64), primary_key=True),
        sa.Column("timestamp", sa.DateTime(timezone=True)),
        sa.Column("actor_id", sa.String(length=128), nullable=False),
        sa.Column("tenant_id", sa.String(length=64), nullable=False),
        sa.Column("action", sa.String(length=64), nullable=False),
        sa.Column("resource_type", sa.String(length=64), nullable=False),
        sa.Column("resource_id", sa.String(length=64), nullable=True),
        sa.Column("request_id", sa.String(length=64), nullable=True),
        sa.Column("source_ip", sa.String(length=64), nullable=True),
        sa.Column("result", sa.String(length=32), nullable=False),
        sa.Column("metadata_json", sa.JSON()),
        sa.Column("prev_hash", sa.String(length=64), nullable=True),
        sa.Column("event_hash", sa.String(length=64), nullable=False),
    )
    op.create_index("ix_audit_events_tenant_id", "audit_events", ["tenant_id"])
    op.create_index("ix_audit_events_timestamp", "audit_events", ["timestamp"])


def downgrade() -> None:
    for table in [
        "audit_events",
        "model_dataset_links",
        "external_cohort_versions",
        "external_cohorts",
        "dataset_cases",
        "dataset_versions",
        "datasets",
        "processing_authorizations",
        "upload_sessions",
        "asset_provenance",
        "genomic_assets",
        "case_status_transitions",
        "genomic_cases",
        "retention_policies",
        "samples",
        "identity_mappings",
        "pseudonymous_patients",
        "tenants",
    ]:
        op.drop_table(table)
