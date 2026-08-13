"""FastAPI dependencies / composition root."""

from __future__ import annotations

from collections.abc import Generator
from dataclasses import dataclass
from typing import Annotated

from fastapi import Depends, Header, Request
from fastapi.security import HTTPAuthorizationCredentials, HTTPBearer
from sqlalchemy.orm import Session

from phelix_vault.application.services.assets import AssetService
from phelix_vault.application.services.audit import AuditService
from phelix_vault.application.services.cases import CaseService
from phelix_vault.application.services.datasets import DatasetService
from phelix_vault.application.services.governance import GovernanceService
from phelix_vault.application.services.identity import IdentityService
from phelix_vault.application.services.retention import RetentionService
from phelix_vault.application.services.samples import SampleService
from phelix_vault.application.services.tenant import TenantService
from phelix_vault.config import Settings, get_settings
from phelix_vault.domain.common import DataResidency
from phelix_vault.infrastructure.kms.base import KeyManagementProvider
from phelix_vault.infrastructure.kms.development import DevelopmentKeyProvider
from phelix_vault.infrastructure.kms.external import ExternalKMSProvider
from phelix_vault.infrastructure.logging import new_request_id
from phelix_vault.infrastructure.security import AuthorizationService, Principal, TokenService
from phelix_vault.infrastructure.storage.base import ObjectStorage
from phelix_vault.infrastructure.storage.local import LocalObjectStorage
from phelix_vault.infrastructure.storage.s3 import S3ObjectStorage

_bearer = HTTPBearer(auto_error=False)


@dataclass
class AppContainer:
    settings: Settings
    engine: object
    session_factory: object
    storage: ObjectStorage
    kms: KeyManagementProvider
    token_service: TokenService
    authz: AuthorizationService


_container: AppContainer | None = None


def set_container(container: AppContainer) -> None:
    global _container
    _container = container


def get_container() -> AppContainer:
    if _container is None:
        raise RuntimeError("Application container not initialized")
    return _container


def get_db() -> Generator[Session, None, None]:
    container = get_container()
    session = container.session_factory()  # type: ignore[operator]
    try:
        yield session
        session.commit()
    except Exception:
        session.rollback()
        raise
    finally:
        session.close()


def get_request_id(
    request: Request,
    x_request_id: Annotated[str | None, Header(alias="X-Request-ID")] = None,
) -> str:
    rid = x_request_id or new_request_id()
    request.state.request_id = rid
    return rid


def get_principal(
    creds: Annotated[HTTPAuthorizationCredentials | None, Depends(_bearer)],
) -> Principal:
    container = get_container()
    if creds is None or not creds.credentials:
        from phelix_vault.domain.errors import AuthenticationError

        raise AuthenticationError("Missing bearer token")
    return container.token_service.verify(creds.credentials)


def build_services(session: Session) -> dict[str, object]:
    container = get_container()
    audit = AuditService(session)
    authz = container.authz
    return {
        "audit": audit,
        "tenants": TenantService(session, audit),
        "identity": IdentityService(session, audit, authz),
        "samples": SampleService(session, audit, authz),
        "cases": CaseService(session, audit, authz),
        "assets": AssetService(
            session,
            audit,
            authz,
            container.storage,
            container.kms,
            container.settings,
        ),
        "datasets": DatasetService(session, audit, authz),
        "governance": GovernanceService(session, audit, authz),
        "retention": RetentionService(session, audit, authz),
    }


def create_container(settings: Settings | None = None) -> AppContainer:
    settings = settings or get_settings()
    from phelix_vault.infrastructure.database import create_db_engine, make_session_factory

    engine = create_db_engine(settings.database_url, echo=settings.debug)
    session_factory = make_session_factory(engine)

    residency = DataResidency(settings.storage_residency)
    storage: ObjectStorage
    if settings.storage_backend == "local":
        storage = LocalObjectStorage(settings.local_storage_root, residency=residency)
    else:
        storage = S3ObjectStorage(
            bucket=settings.s3_bucket,
            region=settings.s3_region,
            residency=residency,
            endpoint_url=settings.s3_endpoint_url,
            access_key=settings.s3_access_key,
            secret_key=settings.s3_secret_key,
        )

    kms: KeyManagementProvider
    if settings.kms_provider == "development":
        kms = DevelopmentKeyProvider(settings.development_kek)
    else:
        kms = ExternalKMSProvider()

    return AppContainer(
        settings=settings,
        engine=engine,
        session_factory=session_factory,
        storage=storage,
        kms=kms,
        token_service=TokenService(settings),
        authz=AuthorizationService(),
    )
