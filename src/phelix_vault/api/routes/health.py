"""Health / readiness — no sensitive information."""

from __future__ import annotations

from fastapi import APIRouter

from phelix_vault import __version__

router = APIRouter(tags=["health"])


@router.get("/health")
def health() -> dict[str, str]:
    return {"status": "ok", "service": "phelix-vault", "version": __version__}


@router.get("/ready")
def ready() -> dict[str, str]:
    # Container presence implies wiring; DB ping left to deployment probes.
    return {"status": "ready"}
