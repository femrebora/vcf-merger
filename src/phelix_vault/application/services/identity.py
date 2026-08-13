"""Identity domain: pseudonymous patients + isolated identity mapping."""

from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import AuditAction, Permission, new_id
from phelix_vault.domain.errors import ConflictError, NotFoundError
from phelix_vault.infrastructure.database.models import IdentityMappingRow, PseudonymousPatientRow
from phelix_vault.infrastructure.security import AuthorizationService, Principal


class IdentityService:
    def __init__(
        self,
        session: Session,
        audit: AuditService,
        authz: AuthorizationService,
    ) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create_patient(
        self,
        principal: Principal,
        *,
        external_subject_id: str | None = None,
        request_id: str | None = None,
    ) -> PseudonymousPatientRow:
        self._authz.require(principal.roles, Permission.PATIENT_CREATE)
        patient = PseudonymousPatientRow(
            id=new_id("PHX-PAT"),
            tenant_id=principal.tenant_id,
        )
        self._session.add(patient)
        self._session.flush()

        if external_subject_id is not None:
            self._authz.require(principal.roles, Permission.IDENTITY_MAPPING_WRITE)
            mapping = IdentityMappingRow(
                id=new_id("PHX-IDM"),
                tenant_id=principal.tenant_id,
                patient_id=patient.id,
                external_subject_id=external_subject_id,
            )
            self._session.add(mapping)

        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.PATIENT_CREATED,
            resource_type="patient",
            resource_id=patient.id,
            request_id=request_id,
            metadata={"has_identity_mapping": external_subject_id is not None},
        )
        self._session.flush()
        return patient

    def get_patient(self, principal: Principal, patient_id: str) -> PseudonymousPatientRow:
        self._authz.require(principal.roles, Permission.PATIENT_READ)
        row = self._session.get(PseudonymousPatientRow, patient_id)
        if not row or row.tenant_id != principal.tenant_id:
            raise NotFoundError("Patient not found")
        return row

    def get_identity_mapping(
        self, principal: Principal, patient_id: str, *, request_id: str | None = None
    ) -> IdentityMappingRow:
        self._authz.require(principal.roles, Permission.IDENTITY_MAPPING_READ)
        row = self._session.scalar(
            select(IdentityMappingRow).where(
                IdentityMappingRow.patient_id == patient_id,
                IdentityMappingRow.tenant_id == principal.tenant_id,
            )
        )
        if not row:
            raise NotFoundError("Identity mapping not found")
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.IDENTITY_MAPPING_ACCESSED,
            resource_type="identity_mapping",
            resource_id=row.id,
            request_id=request_id,
        )
        return row

    def ensure_unique_external(self, tenant_id: str, external_subject_id: str) -> None:
        existing = self._session.scalar(
            select(IdentityMappingRow).where(
                IdentityMappingRow.tenant_id == tenant_id,
                IdentityMappingRow.external_subject_id == external_subject_id,
            )
        )
        if existing:
            raise ConflictError("External subject already mapped")
