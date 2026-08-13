"""Dataset versioning with TRAINING / EXTERNAL_VALIDATION leakage safeguards."""

from __future__ import annotations

from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import (
    AssetType,
    AuditAction,
    DatasetApprovalState,
    DatasetRole,
    DatasetSource,
    MlEligibilityStatus,
    Permission,
    ProcessingPurpose,
    new_id,
    utcnow,
)
from phelix_vault.domain.errors import DatasetLeakageError, NotFoundError, ValidationError
from phelix_vault.infrastructure.database.models import (
    DatasetCaseRow,
    DatasetRow,
    DatasetVersionRow,
    GenomicAssetRow,
    GenomicCaseRow,
)
from phelix_vault.infrastructure.security import AuthorizationService, Principal

# Roles that must never share patients silently.
_CONFLICTING_ROLE_PAIRS: set[frozenset[DatasetRole]] = {
    frozenset({DatasetRole.TRAINING, DatasetRole.EXTERNAL_VALIDATION}),
    frozenset({DatasetRole.TRAINING, DatasetRole.TEST}),
    frozenset({DatasetRole.VALIDATION, DatasetRole.EXTERNAL_VALIDATION}),
}


class DatasetService:
    def __init__(self, session: Session, audit: AuditService, authz: AuthorizationService) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create_dataset(
        self,
        principal: Principal,
        *,
        name: str,
        dataset_role: DatasetRole | str,
        source: DatasetSource | str = DatasetSource.INTERNAL,
        request_id: str | None = None,
    ) -> DatasetRow:
        self._authz.require(principal.roles, Permission.DATASET_MANAGE)
        role = dataset_role.value if isinstance(dataset_role, DatasetRole) else dataset_role
        src = source.value if isinstance(source, DatasetSource) else source
        try:
            DatasetRole(role)
            DatasetSource(src)
        except ValueError as exc:
            raise ValidationError("Invalid dataset role or source") from exc

        row = DatasetRow(
            id=new_id("PHX-DS"),
            tenant_id=principal.tenant_id,
            name=name,
            dataset_role=role,
            source=src,
            created_by=principal.subject,
        )
        self._session.add(row)
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.DATASET_CREATED,
            resource_type="dataset",
            resource_id=row.id,
            request_id=request_id,
            metadata={"role": role, "name": name},
        )
        return row

    def create_version(
        self,
        principal: Principal,
        dataset_id: str,
        *,
        version: str,
        purpose: ProcessingPurpose | str,
        case_ids: list[str],
        excluded_case_ids: list[str] | None = None,
        eligibility_criteria: dict[str, Any] | None = None,
        request_id: str | None = None,
    ) -> DatasetVersionRow:
        self._authz.require(principal.roles, Permission.DATASET_MANAGE)
        dataset = self._session.get(DatasetRow, dataset_id)
        if not dataset or dataset.tenant_id != principal.tenant_id:
            raise NotFoundError("Dataset not found")

        purpose_value = purpose.value if isinstance(purpose, ProcessingPurpose) else purpose
        role = DatasetRole(dataset.dataset_role)

        cases = self._load_cases(principal.tenant_id, case_ids)
        self._assert_cases_eligible(cases, role)
        self._assert_no_patient_leakage(
            tenant_id=principal.tenant_id,
            role=role,
            patient_ids={c.patient_id for c in cases},
            exclude_dataset_id=dataset.id,
        )

        dv = DatasetVersionRow(
            id=new_id("PHX-DSV"),
            tenant_id=principal.tenant_id,
            dataset_id=dataset.id,
            version=version,
            purpose=purpose_value,
            dataset_role=role.value,
            approval_state=DatasetApprovalState.DRAFT.value,
            eligibility_criteria=eligibility_criteria or {},
            excluded_case_ids=excluded_case_ids or [],
            created_by=principal.subject,
        )
        self._session.add(dv)
        self._session.flush()

        for case in cases:
            vcf_id = self._find_vcf_asset(principal.tenant_id, case)
            self._session.add(
                DatasetCaseRow(
                    id=new_id("PHX-DSC"),
                    tenant_id=principal.tenant_id,
                    dataset_version_id=dv.id,
                    case_id=case.id,
                    sample_id=case.sample_id,
                    patient_id=case.patient_id,
                    vcf_asset_id=vcf_id,
                )
            )
        self._session.flush()
        return dv

    def approve_version(
        self,
        principal: Principal,
        dataset_version_id: str,
        *,
        request_id: str | None = None,
    ) -> DatasetVersionRow:
        self._authz.require(principal.roles, Permission.DATASET_APPROVE)
        dv = self._get_version(principal, dataset_version_id)
        if dv.immutable:
            raise ValidationError("Dataset version is immutable")
        # Re-check leakage at approval time.
        patient_ids = {
            row.patient_id
            for row in self._session.scalars(
                select(DatasetCaseRow).where(DatasetCaseRow.dataset_version_id == dv.id)
            )
        }
        self._assert_no_patient_leakage(
            tenant_id=principal.tenant_id,
            role=DatasetRole(dv.dataset_role),
            patient_ids=patient_ids,
            exclude_dataset_id=dv.dataset_id,
        )
        dv.approval_state = DatasetApprovalState.APPROVED.value
        dv.approved_by = principal.subject
        dv.approved_at = utcnow()
        dv.immutable = True
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.DATASET_APPROVED,
            resource_type="dataset_version",
            resource_id=dv.id,
            request_id=request_id,
        )
        return dv

    def export_manifest(
        self,
        principal: Principal,
        dataset_version_id: str,
        *,
        request_id: str | None = None,
    ) -> dict[str, Any]:
        self._authz.require(principal.roles, Permission.DATASET_EXPORT)
        dv = self._get_version(principal, dataset_version_id)
        if dv.approval_state != DatasetApprovalState.APPROVED.value:
            raise ValidationError("Only approved dataset versions may be exported")

        cases = list(
            self._session.scalars(
                select(DatasetCaseRow).where(DatasetCaseRow.dataset_version_id == dv.id)
            )
        )
        manifest = {
            "dataset_id": dv.dataset_id,
            "dataset_version_id": dv.id,
            "version": dv.version,
            "purpose": dv.purpose,
            "dataset_role": dv.dataset_role,
            "created_at": dv.created_at.isoformat() if dv.created_at else None,
            "created_by": dv.created_by,
            "approved_by": dv.approved_by,
            "approved_at": dv.approved_at.isoformat() if dv.approved_at else None,
            "eligibility_criteria": dv.eligibility_criteria,
            "excluded_cases": dv.excluded_case_ids,
            "cases": [
                {
                    "case_id": c.case_id,
                    "sample_id": c.sample_id,
                    "patient_id": c.patient_id,
                    "vcf_asset_id": c.vcf_asset_id,
                    "phenotype_ref": c.phenotype_ref,
                    "truth_label_ref": c.truth_label_ref,
                }
                for c in cases
            ],
        }
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.DATASET_EXPORTED,
            resource_type="dataset_version",
            resource_id=dv.id,
            request_id=request_id,
            metadata={"case_count": len(cases)},
        )
        return manifest

    def patient_overlap_report(
        self, principal: Principal, version_a: str, version_b: str
    ) -> dict[str, Any]:
        self._authz.require(principal.roles, Permission.DATASET_MANAGE)
        a = self._get_version(principal, version_a)
        b = self._get_version(principal, version_b)
        patients_a = self._patients_for_version(a.id)
        patients_b = self._patients_for_version(b.id)
        overlap = sorted(patients_a & patients_b)
        return {
            "version_a": a.id,
            "version_b": b.id,
            "role_a": a.dataset_role,
            "role_b": b.dataset_role,
            "patient_overlap": overlap,
            "overlap_count": len(overlap),
            "conflicting_roles": self._roles_conflict(
                DatasetRole(a.dataset_role), DatasetRole(b.dataset_role)
            ),
        }

    def _get_version(self, principal: Principal, dataset_version_id: str) -> DatasetVersionRow:
        dv = self._session.get(DatasetVersionRow, dataset_version_id)
        if not dv or dv.tenant_id != principal.tenant_id:
            raise NotFoundError("Dataset version not found")
        return dv

    def _load_cases(self, tenant_id: str, case_ids: list[str]) -> list[GenomicCaseRow]:
        if not case_ids:
            raise ValidationError("Dataset version requires at least one case")
        rows = list(
            self._session.scalars(
                select(GenomicCaseRow).where(
                    GenomicCaseRow.tenant_id == tenant_id,
                    GenomicCaseRow.id.in_(case_ids),
                )
            )
        )
        if len(rows) != len(set(case_ids)):
            raise NotFoundError("One or more cases not found in tenant")
        return rows

    def _assert_cases_eligible(self, cases: list[GenomicCaseRow], role: DatasetRole) -> None:
        for case in cases:
            if case.ml_eligibility == MlEligibilityStatus.REVOKED.value:
                raise ValidationError(f"Case {case.id} ML eligibility is REVOKED")
            if role in {
                DatasetRole.TRAINING,
                DatasetRole.VALIDATION,
                DatasetRole.TEST,
                DatasetRole.EXTERNAL_VALIDATION,
            }:
                if case.ml_eligibility != MlEligibilityStatus.APPROVED.value:
                    raise ValidationError(
                        f"Case {case.id} is not APPROVED for ML (status={case.ml_eligibility})"
                    )

    def _assert_no_patient_leakage(
        self,
        *,
        tenant_id: str,
        role: DatasetRole,
        patient_ids: set[str],
        exclude_dataset_id: str,
    ) -> None:
        # Approved versions with conflicting roles that already include these patients.
        versions = list(
            self._session.scalars(
                select(DatasetVersionRow).where(
                    DatasetVersionRow.tenant_id == tenant_id,
                    DatasetVersionRow.approval_state == DatasetApprovalState.APPROVED.value,
                    DatasetVersionRow.dataset_id != exclude_dataset_id,
                )
            )
        )
        for dv in versions:
            other_role = DatasetRole(dv.dataset_role)
            if not self._roles_conflict(role, other_role):
                continue
            existing = self._patients_for_version(dv.id)
            overlap = patient_ids & existing
            if overlap:
                raise DatasetLeakageError(
                    f"Patient overlap between {role.value} and {other_role.value}: "
                    f"{len(overlap)} patient(s)"
                )

    def _patients_for_version(self, dataset_version_id: str) -> set[str]:
        return {
            row.patient_id
            for row in self._session.scalars(
                select(DatasetCaseRow).where(
                    DatasetCaseRow.dataset_version_id == dataset_version_id
                )
            )
        }

    @staticmethod
    def _roles_conflict(a: DatasetRole, b: DatasetRole) -> bool:
        if a == b:
            return False
        return frozenset({a, b}) in _CONFLICTING_ROLE_PAIRS

    def _find_vcf_asset(self, tenant_id: str, case: GenomicCaseRow) -> str | None:
        asset = self._session.scalar(
            select(GenomicAssetRow).where(
                GenomicAssetRow.tenant_id == tenant_id,
                GenomicAssetRow.case_id == case.id,
                GenomicAssetRow.asset_type == AssetType.VCF.value,
            )
        )
        if asset:
            return asset.id
        asset = self._session.scalar(
            select(GenomicAssetRow).where(
                GenomicAssetRow.tenant_id == tenant_id,
                GenomicAssetRow.sample_id == case.sample_id,
                GenomicAssetRow.asset_type == AssetType.VCF.value,
            )
        )
        return asset.id if asset else None
