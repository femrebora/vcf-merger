"""Runtime configuration via environment variables."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Literal

from pydantic import Field
from pydantic_settings import BaseSettings, SettingsConfigDict


class Settings(BaseSettings):
    model_config = SettingsConfigDict(env_file=".env", env_file_encoding="utf-8", extra="ignore")

    app_name: str = "pHelix Vault"
    environment: Literal["development", "test", "production"] = "development"
    debug: bool = False

    api_prefix: str = "/api/v1"
    cors_allowed_origins: str = "http://localhost:3000"

    database_url: str = "sqlite+pysqlite:////tmp/phelix-vault.db"

    jwt_secret: str = "CHANGE_ME_DEV_ONLY_NOT_FOR_PRODUCTION"
    jwt_algorithm: str = "HS256"
    jwt_audience: str = "phelix-vault"
    jwt_issuer: str = "phelix-dev"

    storage_backend: Literal["local", "s3"] = "local"
    local_storage_root: Path = Path("/tmp/phelix-vault-objects")
    s3_endpoint_url: str | None = None
    s3_access_key: str | None = None
    s3_secret_key: str | None = None
    s3_bucket: str = "phelix-vault"
    s3_region: str = "us-east-1"
    storage_residency: Literal["TR", "EU", "US", "OTHER"] = "TR"

    upload_session_ttl_seconds: int = 3600
    signed_url_ttl_seconds: int = 300
    max_upload_bytes: int = 50 * 1024 * 1024 * 1024  # 50 GiB metadata limit

    kms_provider: Literal["development", "external"] = "development"
    development_kek: str = "dev-only-kek-not-for-production"

    enable_rls: bool = False
    request_id_header: str = "X-Request-ID"

    rate_limit_per_minute: int = Field(default=600, ge=1)


@lru_cache
def get_settings() -> Settings:
    return Settings()
