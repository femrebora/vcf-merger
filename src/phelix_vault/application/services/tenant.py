"""Tenant application service."""

from __future__ import annotations

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import AuditAction, DataResidency, new_id, utcnow
from phelix_vault.domain.errors import ConflictError, NotFoundError
from phelix_vault.infrastructure.database.models import TenantRow


class TenantService:
    def __init__(self, session: Session, audit: AuditService) -> None:
        self._session = session
        self._audit = audit

    def create(
        self,
        *,
        name: str,
        required_residency: DataResidency | str = DataResidency.TR,
        actor_id: str,
        request_id: str | None = None,
    ) -> TenantRow:
        existing = self._session.scalar(select(TenantRow).where(TenantRow.name == name))
        if existing:
            raise ConflictError("Tenant name already exists")
        residency = (
            required_residency.value
            if isinstance(required_residency, DataResidency)
            else required_residency
        )
        row = TenantRow(id=new_id("PHX-TNT"), name=name, required_residency=residency)
        self._session.add(row)
        self._session.flush()
        self._audit.record(
            actor_id=actor_id,
            tenant_id=row.id,
            action=AuditAction.TENANT_CREATED,
            resource_type="tenant",
            resource_id=row.id,
            request_id=request_id,
            metadata={"name": name},
        )
        return row

    def get(self, tenant_id: str) -> TenantRow:
        row = self._session.get(TenantRow, tenant_id)
        if not row:
            raise NotFoundError("Tenant not found")
        return row

    def touch(self, tenant_id: str) -> TenantRow:
        row = self.get(tenant_id)
        row.updated_at = utcnow()
        row.version += 1
        self._session.flush()
        return row
