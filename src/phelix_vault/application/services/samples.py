"""Sample registry service."""

from __future__ import annotations

from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import AuditAction, Permission, new_id
from phelix_vault.domain.errors import NotFoundError
from phelix_vault.infrastructure.database.models import PseudonymousPatientRow, SampleRow
from phelix_vault.infrastructure.security import AuthorizationService, Principal


class SampleService:
    def __init__(self, session: Session, audit: AuditService, authz: AuthorizationService) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create(
        self,
        principal: Principal,
        *,
        patient_id: str,
        sample_type: str = "WES",
        collection_date: str | None = None,
        request_id: str | None = None,
    ) -> SampleRow:
        self._authz.require(principal.roles, Permission.SAMPLE_MANAGE)
        patient = self._session.get(PseudonymousPatientRow, patient_id)
        if not patient or patient.tenant_id != principal.tenant_id:
            raise NotFoundError("Patient not found")
        row = SampleRow(
            id=new_id("PHX-SMP"),
            tenant_id=principal.tenant_id,
            patient_id=patient_id,
            sample_type=sample_type,
            collection_date=collection_date,
        )
        self._session.add(row)
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.SAMPLE_CREATED,
            resource_type="sample",
            resource_id=row.id,
            request_id=request_id,
        )
        return row

    def get(self, principal: Principal, sample_id: str) -> SampleRow:
        self._authz.require(principal.roles, Permission.SAMPLE_MANAGE)
        row = self._session.get(SampleRow, sample_id)
        if not row or row.tenant_id != principal.tenant_id:
            raise NotFoundError("Sample not found")
        return row
