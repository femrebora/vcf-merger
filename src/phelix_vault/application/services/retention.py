"""Retention policy evaluation (no automatic deletion in development)."""

from __future__ import annotations

from datetime import timedelta

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.application.services.audit import AuditService
from phelix_vault.domain.common import (
    AssetStatus,
    AuditAction,
    Permission,
    RetentionAction,
    as_utc,
    new_id,
    utcnow,
)
from phelix_vault.domain.errors import ValidationError
from phelix_vault.infrastructure.database.models import GenomicAssetRow, RetentionPolicyRow
from phelix_vault.infrastructure.security import AuthorizationService, Principal


class RetentionService:
    def __init__(self, session: Session, audit: AuditService, authz: AuthorizationService) -> None:
        self._session = session
        self._audit = audit
        self._authz = authz

    def create_policy(
        self,
        principal: Principal,
        *,
        name: str,
        asset_type: str,
        retention_period_days: int,
        action: RetentionAction | str,
        legal_hold_supported: bool = True,
        request_id: str | None = None,
    ) -> RetentionPolicyRow:
        self._authz.require(principal.roles, Permission.RETENTION_MANAGE)
        if retention_period_days < 1:
            raise ValidationError("retention_period_days must be >= 1")
        action_value = action.value if isinstance(action, RetentionAction) else action
        row = RetentionPolicyRow(
            id=new_id("PHX-RET"),
            tenant_id=principal.tenant_id,
            name=name,
            asset_type=asset_type,
            retention_period_days=retention_period_days,
            action=action_value,
            legal_hold_supported=legal_hold_supported,
        )
        self._session.add(row)
        self._session.flush()
        self._audit.record(
            actor_id=principal.subject,
            tenant_id=principal.tenant_id,
            action=AuditAction.RETENTION_POLICY_CHANGED,
            resource_type="retention_policy",
            resource_id=row.id,
            request_id=request_id,
            metadata={"name": name, "action": action_value},
        )
        return row

    def evaluate(self, principal: Principal) -> list[dict[str, str | None]]:
        """Return assets that would be due for retention action — does not delete."""
        self._authz.require(principal.roles, Permission.RETENTION_MANAGE)
        policies = list(
            self._session.scalars(
                select(RetentionPolicyRow).where(
                    RetentionPolicyRow.tenant_id == principal.tenant_id
                )
            )
        )
        due: list[dict[str, str | None]] = []
        now = utcnow()
        for policy in policies:
            assets = self._session.scalars(
                select(GenomicAssetRow).where(
                    GenomicAssetRow.tenant_id == principal.tenant_id,
                    GenomicAssetRow.asset_type == policy.asset_type,
                    GenomicAssetRow.status == AssetStatus.AVAILABLE.value,
                )
            )
            for asset in assets:
                ref = asset.ingested_at or asset.created_at
                if ref and as_utc(ref) + timedelta(days=policy.retention_period_days) <= now:
                    due.append(
                        {
                            "asset_id": asset.id,
                            "policy_id": policy.id,
                            "action": policy.action,
                            "asset_type": asset.asset_type,
                        }
                    )
        return due
