"""Centralized RBAC policy layer."""

from __future__ import annotations

from phelix_vault.domain.common import Permission, Role
from phelix_vault.domain.errors import AuthorizationError

ROLE_PERMISSIONS: dict[Role, frozenset[Permission]] = {
    Role.SYSTEM_ADMIN: frozenset(Permission),
    Role.TENANT_ADMIN: frozenset(
        {
            Permission.TENANT_MANAGE,
            Permission.PATIENT_CREATE,
            Permission.PATIENT_READ,
            Permission.IDENTITY_MAPPING_READ,
            Permission.IDENTITY_MAPPING_WRITE,
            Permission.SAMPLE_MANAGE,
            Permission.ASSET_UPLOAD,
            Permission.ASSET_READ_META,
            Permission.ASSET_DOWNLOAD,
            Permission.ASSET_DELETE,
            Permission.CASE_CURATE,
            Permission.ML_ELIGIBILITY_MANAGE,
            Permission.DATASET_MANAGE,
            Permission.DATASET_APPROVE,
            Permission.DATASET_EXPORT,
            Permission.GOVERNANCE_MANAGE,
            Permission.RETENTION_MANAGE,
            Permission.AUDIT_READ,
        }
    ),
    Role.CLINICAL_GENETICIST: frozenset(
        {
            Permission.PATIENT_READ,
            Permission.SAMPLE_MANAGE,
            Permission.ASSET_READ_META,
            Permission.ASSET_DOWNLOAD,
            Permission.CASE_CURATE,
            Permission.ML_ELIGIBILITY_MANAGE,
        }
    ),
    Role.BIOINFORMATICIAN: frozenset(
        {
            Permission.SAMPLE_MANAGE,
            Permission.ASSET_UPLOAD,
            Permission.ASSET_READ_META,
            Permission.ASSET_DOWNLOAD,
            Permission.CASE_CURATE,
        }
    ),
    Role.LAB_TECHNICIAN: frozenset(
        {
            Permission.SAMPLE_MANAGE,
            Permission.ASSET_UPLOAD,
            Permission.ASSET_READ_META,
        }
    ),
    Role.DATA_CURATOR: frozenset(
        {
            Permission.PATIENT_READ,
            Permission.SAMPLE_MANAGE,
            Permission.ASSET_READ_META,
            Permission.CASE_CURATE,
            Permission.ML_ELIGIBILITY_MANAGE,
            Permission.DATASET_MANAGE,
        }
    ),
    Role.ML_RESEARCHER: frozenset(
        {
            Permission.ASSET_READ_META,
            Permission.DATASET_MANAGE,
            Permission.DATASET_EXPORT,
        }
    ),
    Role.AUDITOR: frozenset(
        {
            Permission.AUDIT_READ,
            Permission.ASSET_READ_META,
            Permission.PATIENT_READ,
        }
    ),
}


class AuthorizationService:
    """Evaluate permissions; never scatter role checks in controllers."""

    def permissions_for(self, roles: list[Role] | list[str]) -> set[Permission]:
        perms: set[Permission] = set()
        for role in roles:
            r = Role(role) if isinstance(role, str) else role
            perms |= set(ROLE_PERMISSIONS.get(r, frozenset()))
        return perms

    def has_permission(self, roles: list[Role] | list[str], permission: Permission) -> bool:
        return permission in self.permissions_for(roles)

    def require(
        self,
        roles: list[Role] | list[str],
        permission: Permission,
        *,
        detail: str = "Insufficient permissions",
    ) -> None:
        if not self.has_permission(roles, permission):
            raise AuthorizationError(detail)

    def require_tenant(self, actor_tenant_id: str, resource_tenant_id: str) -> None:
        if actor_tenant_id != resource_tenant_id:
            raise AuthorizationError("Cross-tenant access denied")
