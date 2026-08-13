"""Append-oriented audit logging (tamper-evident hash chain ready)."""

from __future__ import annotations

import hashlib
import json
from typing import Any

from sqlalchemy import select
from sqlalchemy.orm import Session

from phelix_vault.domain.common import AuditAction, new_id, utcnow
from phelix_vault.infrastructure.database.models import AuditEventRow


class AuditService:
    def __init__(self, session: Session) -> None:
        self._session = session

    def record(
        self,
        *,
        actor_id: str,
        tenant_id: str,
        action: AuditAction | str,
        resource_type: str,
        resource_id: str | None = None,
        request_id: str | None = None,
        source_ip: str | None = None,
        result: str = "SUCCESS",
        metadata: dict[str, Any] | None = None,
    ) -> AuditEventRow:
        action_value = action.value if isinstance(action, AuditAction) else action
        meta = metadata or {}
        # Never store sensitive payloads in audit metadata.
        for banned in ("vcf", "fastq", "sequence", "patient_name", "tc_no", "report_text"):
            meta.pop(banned, None)

        prev = self._session.scalar(
            select(AuditEventRow)
            .where(AuditEventRow.tenant_id == tenant_id)
            .order_by(AuditEventRow.timestamp.desc())
            .limit(1)
        )
        prev_hash = prev.event_hash if prev else None
        event_id = new_id("PHX-AUD")
        ts = utcnow()
        material = json.dumps(
            {
                "id": event_id,
                "timestamp": ts.isoformat(),
                "actor_id": actor_id,
                "tenant_id": tenant_id,
                "action": action_value,
                "resource_type": resource_type,
                "resource_id": resource_id,
                "result": result,
                "prev_hash": prev_hash,
            },
            sort_keys=True,
        )
        event_hash = hashlib.sha256(material.encode("utf-8")).hexdigest()
        row = AuditEventRow(
            id=event_id,
            timestamp=ts,
            actor_id=actor_id,
            tenant_id=tenant_id,
            action=action_value,
            resource_type=resource_type,
            resource_id=resource_id,
            request_id=request_id,
            source_ip=source_ip,
            result=result,
            metadata_json=meta,
            prev_hash=prev_hash,
            event_hash=event_hash,
        )
        self._session.add(row)
        self._session.flush()
        return row

    def list_for_tenant(self, tenant_id: str, *, limit: int = 100) -> list[AuditEventRow]:
        return list(
            self._session.scalars(
                select(AuditEventRow)
                .where(AuditEventRow.tenant_id == tenant_id)
                .order_by(AuditEventRow.timestamp.desc())
                .limit(limit)
            )
        )
