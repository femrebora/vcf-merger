"""Domain-specific errors (never leak secrets/paths to clients)."""

from __future__ import annotations


class DomainError(Exception):
    code: str = "domain_error"
    http_status: int = 400

    def __init__(self, message: str = "Domain error") -> None:
        self.message = message
        super().__init__(message)


class NotFoundError(DomainError):
    code = "not_found"
    http_status = 404


class ConflictError(DomainError):
    code = "conflict"
    http_status = 409


class AuthorizationError(DomainError):
    code = "forbidden"
    http_status = 403


class AuthenticationError(DomainError):
    code = "unauthorized"
    http_status = 401


class ValidationError(DomainError):
    code = "validation_error"
    http_status = 422


class DatasetLeakageError(ConflictError):
    code = "dataset_leakage"


class IntegrityError(DomainError):
    code = "integrity_error"
    http_status = 400


class ResidencyViolationError(DomainError):
    code = "residency_violation"
    http_status = 403
