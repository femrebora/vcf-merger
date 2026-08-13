"""Unit tests for RBAC policy layer."""

from __future__ import annotations

import pytest

from phelix_vault.domain.common import Permission, Role
from phelix_vault.domain.errors import AuthorizationError
from phelix_vault.infrastructure.security.rbac import AuthorizationService


def test_ml_researcher_cannot_access_identity_mapping() -> None:
    authz = AuthorizationService()
    assert not authz.has_permission([Role.ML_RESEARCHER], Permission.IDENTITY_MAPPING_READ)
    with pytest.raises(AuthorizationError):
        authz.require([Role.ML_RESEARCHER], Permission.IDENTITY_MAPPING_READ)


def test_auditor_cannot_download_assets() -> None:
    authz = AuthorizationService()
    assert not authz.has_permission([Role.AUDITOR], Permission.ASSET_DOWNLOAD)


def test_lab_technician_cannot_approve_datasets() -> None:
    authz = AuthorizationService()
    assert not authz.has_permission([Role.LAB_TECHNICIAN], Permission.DATASET_APPROVE)


def test_tenant_admin_has_core_permissions() -> None:
    authz = AuthorizationService()
    for perm in (
        Permission.ASSET_UPLOAD,
        Permission.ML_ELIGIBILITY_MANAGE,
        Permission.DATASET_APPROVE,
        Permission.AUDIT_READ,
    ):
        assert authz.has_permission([Role.TENANT_ADMIN], perm)


def test_cross_tenant_require() -> None:
    authz = AuthorizationService()
    with pytest.raises(AuthorizationError):
        authz.require_tenant("t1", "t2")
