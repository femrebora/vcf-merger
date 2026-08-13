"""Governance metadata — records decisions; does not assert legal compliance."""

from __future__ import annotations

from datetime import datetime

from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import Permission, ProcessingPurpose, new_id, utcnow
from phelix_vault.domain.errors import NotFoundError, ValidationError
from phelix_vault.infrastructure.database.models import GenomicCaseRow, ProcessingAuthorizationRow
from phelix_vault.infrastructure.security import AuthorizationService, Principal


class GovernanceService:
    def __init__(self, session: Session, audit: AuditService, authz: AuthorizationService) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create_authorization(
        self,
        principal: Principal,
        *,
        case_id: str,
        processing_purpose: ProcessingPurpose | str,
        legal_basis_code: str,
        authorization_source: str,
        effective_from: datetime | None = None,
        effective_until: datetime | None = None,
        restrictions: str | None = None,
        data_controller: str | None = None,
        data_processor: str | None = None,
        international_transfer_allowed: bool = False,
        ml_training_allowed: bool = False,
        research_allowed: bool = False,
        reanalysis_allowed: bool = True,
        notes: str | None = None,
    ) -> ProcessingAuthorizationRow:
        self._authz.require(principal.roles, Permission.GOVERNANCE_MANAGE)
        case = self._session.get(GenomicCaseRow, case_id)
        if not case or case.tenant_id != principal.tenant_id:
            raise NotFoundError("Case not found")

        purpose = (
            processing_purpose.value
            if isinstance(processing_purpose, ProcessingPurpose)
            else processing_purpose
        )
        try:
            ProcessingPurpose(purpose)
        except ValueError as exc:
            raise ValidationError("Invalid processing purpose") from exc

        if purpose == ProcessingPurpose.ML_TRAINING.value and not ml_training_allowed:
            raise ValidationError("ML_TRAINING purpose requires ml_training_allowed=true")

        row = ProcessingAuthorizationRow(
            id=new_id("PHX-AUTH"),
            tenant_id=principal.tenant_id,
            case_id=case_id,
            processing_purpose=purpose,
            legal_basis_code=legal_basis_code,
            authorization_source=authorization_source,
            effective_from=effective_from or utcnow(),
            effective_until=effective_until,
            restrictions=restrictions,
            data_controller=data_controller,
            data_processor=data_processor,
            international_transfer_allowed=international_transfer_allowed,
            ml_training_allowed=ml_training_allowed,
            research_allowed=research_allowed,
            reanalysis_allowed=reanalysis_allowed,
            approved_by=principal.subject,
            notes=notes,
        )
        self._session.add(row)
        self._session.flush()
        return row
