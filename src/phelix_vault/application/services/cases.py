"""Genomic case curation and ML eligibility — explicit transitions only."""

from __future__ import annotations

from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import (
    AuditAction,
    CurationStatus,
    DataZone,
    MlEligibilityStatus,
    Permission,
    ProcessingPurpose,
    new_id,
    utcnow,
)
from phelix_vault.domain.errors import NotFoundError, ValidationError
from phelix_vault.infrastructure.database.models import (
    CaseStatusTransitionRow,
    GenomicCaseRow,
    PseudonymousPatientRow,
    SampleRow,
)
from phelix_vault.infrastructure.security import AuthorizationService, Principal

_CURATION_ORDER = [
    CurationStatus.RAW,
    CurationStatus.TECHNICALLY_VALID,
    CurationStatus.INTERPRETATION_UNKNOWN,
    CurationStatus.PARTIALLY_CURATED,
    CurationStatus.EXPERT_REVIEWED,
    CurationStatus.GOLD_STANDARD,
]


class CaseService:
    def __init__(self, session: Session, audit: AuditService, authz: AuthorizationService) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create(
        self,
        principal: Principal,
        *,
        patient_id: str,
        sample_id: str,
        purposes: list[ProcessingPurpose | str] | None = None,
        request_id: str | None = None,
    ) -> GenomicCaseRow:
        self._authz.require(principal.roles, Permission.CASE_CURATE)
        patient = self._session.get(PseudonymousPatientRow, patient_id)
        sample = self._session.get(SampleRow, sample_id)
        if (
            not patient
            or not sample
            or patient.tenant_id != principal.tenant_id
            or sample.tenant_id != principal.tenant_id
            or sample.patient_id != patient_id
        ):
            raise NotFoundError("Patient or sample not found")

        purpose_values = [
            p.value if isinstance(p, ProcessingPurpose) else p for p in (purposes or [])
        ]
        # Clinical storage purpose may be present; ML is NEVER implied.
        if ProcessingPurpose.ML_TRAINING.value in purpose_values:
            raise ValidationError(
                "ML_TRAINING cannot be set at case creation; use explicit ML eligibility approval"
            )

        row = GenomicCaseRow(
            id=new_id("PHX-CASE"),
            tenant_id=principal.tenant_id,
            patient_id=patient_id,
            sample_id=sample_id,
            curation_status=CurationStatus.RAW.value,
            ml_eligibility=MlEligibilityStatus.NOT_EVALUATED.value,
            data_zone=DataZone.RAW_CLINICAL.value,
            processing_purposes=purpose_values or [ProcessingPurpose.CLINICAL_DIAGNOSIS.value],
        )
        self._session.add(row)
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.CASE_CREATED,
            resource_type="case",
            resource_id=row.id,
            request_id=request_id,
        )
        return row

    def get(self, principal: Principal, case_id: str) -> GenomicCaseRow:
        row = self._session.get(GenomicCaseRow, case_id)
        if not row or row.tenant_id != principal.tenant_id:
            raise NotFoundError("Case not found")
        return row

    def _transition(
        self,
        principal: Principal,
        case: GenomicCaseRow,
        *,
        field_name: str,
        new_state: str,
        reason: str | None,
        request_id: str | None,
        audit_action: AuditAction,
    ) -> GenomicCaseRow:
        previous = getattr(case, field_name)
        if previous == new_state:
            return case
        setattr(case, field_name, new_state)
        case.updated_at = utcnow()
        case.version += 1
        self._session.add(
            CaseStatusTransitionRow(
                id=new_id("PHX-TRN"),
                tenant_id=principal.tenant_id,
                case_id=case.id,
                field_name=field_name,
                previous_state=previous,
                new_state=new_state,
                reason=reason,
                actor_id=principal.subject,
            )
        )
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=audit_action,
            resource_type="case",
            resource_id=case.id,
            request_id=request_id,
            metadata={
                "field": field_name,
                "previous": previous,
                "new": new_state,
                "reason": reason,
            },
        )
        self._session.flush()
        return case

    def set_curation_status(
        self,
        principal: Principal,
        case_id: str,
        *,
        status: CurationStatus | str,
        reason: str | None = None,
        request_id: str | None = None,
    ) -> GenomicCaseRow:
        self._authz.require(principal.roles, Permission.CASE_CURATE)
        case = self.get(principal, case_id)
        new_status = status.value if isinstance(status, CurationStatus) else status
        try:
            CurationStatus(new_status)
        except ValueError as exc:
            raise ValidationError("Invalid curation status") from exc

        # Explicit promotion — do not silently skip, but allow documented non-linear sets
        # with reason for expert override.
        if reason is None:
            current = CurationStatus(case.curation_status)
            target = CurationStatus(new_status)
            if current in _CURATION_ORDER and target in _CURATION_ORDER:
                if _CURATION_ORDER.index(target) > _CURATION_ORDER.index(current) + 1:
                    raise ValidationError(
                        "Non-sequential curation promotion requires an explicit reason"
                    )

        case = self._transition(
            principal,
            case,
            field_name="curation_status",
            new_state=new_status,
            reason=reason,
            request_id=request_id,
            audit_action=AuditAction.CURATION_STATUS_CHANGED,
        )
        if new_status in {
            CurationStatus.PARTIALLY_CURATED.value,
            CurationStatus.EXPERT_REVIEWED.value,
            CurationStatus.GOLD_STANDARD.value,
        }:
            case.data_zone = DataZone.CURATED.value
        elif new_status == CurationStatus.TECHNICALLY_VALID.value:
            case.data_zone = DataZone.PROCESSED_CLINICAL.value
        self._session.flush()
        return case

    def set_ml_eligibility(
        self,
        principal: Principal,
        case_id: str,
        *,
        status: MlEligibilityStatus | str,
        reason: str | None = None,
        request_id: str | None = None,
    ) -> GenomicCaseRow:
        self._authz.require(principal.roles, Permission.ML_ELIGIBILITY_MANAGE)
        case = self.get(principal, case_id)
        new_status = status.value if isinstance(status, MlEligibilityStatus) else status
        try:
            MlEligibilityStatus(new_status)
        except ValueError as exc:
            raise ValidationError("Invalid ML eligibility status") from exc

        if new_status == MlEligibilityStatus.APPROVED.value:
            if case.ml_eligibility == MlEligibilityStatus.REVOKED.value and not reason:
                raise ValidationError("Re-approval after revocation requires a reason")
            # Approving for ML does not assign dataset role; that is a separate governance step.
            purposes = list(case.processing_purposes or [])
            if ProcessingPurpose.ML_TRAINING.value not in purposes:
                purposes.append(ProcessingPurpose.ML_TRAINING.value)
                case.processing_purposes = purposes
            case.data_zone = DataZone.ML_APPROVED.value

        return self._transition(
            principal,
            case,
            field_name="ml_eligibility",
            new_state=new_status,
            reason=reason,
            request_id=request_id,
            audit_action=AuditAction.ML_ELIGIBILITY_CHANGED,
        )

    def eligibility_report(self, principal: Principal, case_id: str) -> dict[str, object]:
        case = self.get(principal, case_id)
        ml_eligible = case.ml_eligibility == MlEligibilityStatus.APPROVED.value
        pipeline_validation_eligible = bool(
            case.technical_completeness and case.pipeline_provenance_known
        )
        return {
            "case_id": case.id,
            "fastq_available": case.technical_completeness,  # populated by callers/inventory
            "vcf_available": case.technical_completeness,
            "report_available": case.report_available,
            "hpo_available": case.phenotype_completeness,
            "confirmed_diagnosis": case.diagnosis_available,
            "expert_reviewed": case.expert_review_count > 0,
            "ml_training": "ELIGIBLE" if ml_eligible else "NOT ELIGIBLE",
            "pipeline_validation": "ELIGIBLE" if pipeline_validation_eligible else "NOT ELIGIBLE",
            "ml_eligibility_status": case.ml_eligibility,
            "curation_status": case.curation_status,
        }
