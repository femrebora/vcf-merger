"""Shared pytest fixtures for pHelix Vault."""

from __future__ import annotations

import os
from collections.abc import Generator
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

# Configure test environment before app import side effects.
os.environ.setdefault("ENVIRONMENT", "test")
os.environ.setdefault("DEBUG", "false")
os.environ.setdefault("JWT_SECRET", "test-jwt-secret-not-for-production")
os.environ.setdefault("STORAGE_BACKEND", "local")
os.environ.setdefault("KMS_PROVIDER", "development")
os.environ.setdefault("DEVELOPMENT_KEK", "test-kek")
os.environ.setdefault("STORAGE_RESIDENCY", "TR")
os.environ.setdefault("CORS_ALLOWED_ORIGINS", "http://localhost:3000")

SYNTHETIC = Path(__file__).resolve().parent / "fixtures" / "synthetic"


@pytest.fixture
def synthetic_dir() -> Path:
    return SYNTHETIC


@pytest.fixture
def client(tmp_path: Path) -> Generator[TestClient, None, None]:
    os.environ["ENVIRONMENT"] = "test"
    os.environ["DATABASE_URL"] = f"sqlite+pysqlite:///{tmp_path / 'vault-test.db'}"
    os.environ["LOCAL_STORAGE_ROOT"] = str(tmp_path / "objects")
    os.environ["JWT_SECRET"] = "test-jwt-secret-not-for-production"
    os.environ["STORAGE_BACKEND"] = "local"
    os.environ["STORAGE_RESIDENCY"] = "TR"
    os.environ["KMS_PROVIDER"] = "development"

    from phelix_vault.config.settings import get_settings

    get_settings.cache_clear()

    from phelix_vault.main import create_app

    app = create_app()
    with TestClient(app) as test_client:
        yield test_client

    get_settings.cache_clear()
