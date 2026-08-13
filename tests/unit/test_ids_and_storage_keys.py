"""Pseudonym and object-key safety."""

from __future__ import annotations

import pytest

from phelix_vault.domain.common import new_id
from phelix_vault.domain.errors import ValidationError
from phelix_vault.infrastructure.storage.base import validate_object_key


def test_opaque_prefixed_ids() -> None:
    pid = new_id("PHX-PAT")
    assert pid.startswith("PHX-PAT-")
    assert " " not in pid
    assert len(pid) > 16


def test_reject_path_traversal_keys() -> None:
    with pytest.raises(ValidationError):
        validate_object_key("../etc/passwd")
    with pytest.raises(ValidationError):
        validate_object_key("/absolute")
    with pytest.raises(ValidationError):
        validate_object_key("a\\b")


def test_accept_tenant_scoped_key() -> None:
    key = validate_object_key("tenants/PHX-TNT-1/assets/PHX-AST-1/vcf")
    assert key.startswith("tenants/")
