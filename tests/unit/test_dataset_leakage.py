"""Dataset leakage invariants at the domain service level."""

from __future__ import annotations

from phelix_vault.application.services.datasets import DatasetService
from phelix_vault.domain.common import DatasetRole


def test_conflicting_role_pairs() -> None:
    assert DatasetService._roles_conflict(DatasetRole.TRAINING, DatasetRole.EXTERNAL_VALIDATION)
    assert DatasetService._roles_conflict(DatasetRole.TRAINING, DatasetRole.TEST)
    assert not DatasetService._roles_conflict(DatasetRole.TRAINING, DatasetRole.TRAINING)
    assert not DatasetService._roles_conflict(DatasetRole.BENCHMARK, DatasetRole.TRAINING)
