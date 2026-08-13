"""FastAPI application entrypoint."""

from __future__ import annotations

from collections.abc import AsyncIterator
from contextlib import asynccontextmanager

from fastapi import FastAPI, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import JSONResponse

from phelix_vault import __version__
from phelix_vault.api.dependencies import create_container, set_container
from phelix_vault.api.routes import (
    assets,
    audit,
    auth_dev,
    cases,
    datasets,
    governance,
    health,
    patients,
    retention,
    samples,
    tenants,
)
from phelix_vault.config import get_settings
from phelix_vault.domain.errors import DomainError
from phelix_vault.infrastructure.database import Base
from phelix_vault.infrastructure.logging import configure_logging, get_logger

logger = get_logger(__name__)


@asynccontextmanager
async def lifespan(app: FastAPI) -> AsyncIterator[None]:
    settings = get_settings()
    configure_logging("DEBUG" if settings.debug else "INFO")
    container = create_container(settings)
    set_container(container)
    # Create schema for MVP (Alembic migrations also provided).
    Base.metadata.create_all(bind=container.engine)  # type: ignore[arg-type]
    app.state.container = container
    logger.info("phelix-vault started", extra={"request_id": "-"})
    yield


def create_app() -> FastAPI:
    settings = get_settings()
    app = FastAPI(
        title=settings.app_name,
        version=__version__,
        description=(
            "Secure genomic data governance and storage for pHelix. "
            "Upload ≠ permission to train. Clinical storage ≠ ML training dataset."
        ),
        lifespan=lifespan,
    )

    origins = [o.strip() for o in settings.cors_allowed_origins.split(",") if o.strip()]
    app.add_middleware(
        CORSMiddleware,
        allow_origins=origins or ["http://localhost:3000"],
        allow_credentials=True,
        allow_methods=["GET", "POST", "PUT", "PATCH", "DELETE"],
        allow_headers=["Authorization", "Content-Type", "X-Request-ID"],
    )

    @app.middleware("http")
    async def security_headers(request: Request, call_next):  # type: ignore[no-untyped-def]
        response = await call_next(request)
        response.headers["X-Content-Type-Options"] = "nosniff"
        response.headers["X-Frame-Options"] = "DENY"
        response.headers["Referrer-Policy"] = "no-referrer"
        response.headers["Cache-Control"] = "no-store"
        rid = getattr(request.state, "request_id", None)
        if rid:
            response.headers["X-Request-ID"] = rid
        return response

    @app.exception_handler(DomainError)
    async def domain_error_handler(request: Request, exc: DomainError) -> JSONResponse:
        rid = getattr(request.state, "request_id", None)
        return JSONResponse(
            status_code=exc.http_status,
            content={"error": exc.code, "message": exc.message, "request_id": rid},
        )

    prefix = settings.api_prefix
    app.include_router(health.router)
    app.include_router(auth_dev.router, prefix=prefix)
    app.include_router(tenants.router, prefix=prefix)
    app.include_router(patients.router, prefix=prefix)
    app.include_router(samples.router, prefix=prefix)
    app.include_router(cases.router, prefix=prefix)
    app.include_router(assets.router, prefix=prefix)
    app.include_router(datasets.router, prefix=prefix)
    app.include_router(governance.router, prefix=prefix)
    app.include_router(retention.router, prefix=prefix)
    app.include_router(audit.router, prefix=prefix)
    return app


app = create_app()


def run() -> None:
    import uvicorn

    settings = get_settings()
    uvicorn.run(
        "phelix_vault.main:app",
        host="0.0.0.0",
        port=8000,
        reload=settings.environment == "development",
    )


if __name__ == "__main__":
    run()
