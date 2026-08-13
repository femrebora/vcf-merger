"""JWT authentication (OIDC-ready; no password store)."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import UTC, datetime, timedelta
from typing import Any

from jose import JWTError, jwt

from phelix_vault.config import Settings
from phelix_vault.domain.common import Role
from phelix_vault.domain.errors import AuthenticationError


@dataclass(slots=True)
class Principal:
    subject: str
    tenant_id: str
    roles: list[Role]
    email: str | None = None


class TokenService:
    def __init__(self, settings: Settings) -> None:
        self._settings = settings

    def create_access_token(
        self,
        *,
        subject: str,
        tenant_id: str,
        roles: list[Role | str],
        email: str | None = None,
        expires_minutes: int = 60,
    ) -> str:
        now = datetime.now(UTC)
        role_values = [r.value if isinstance(r, Role) else r for r in roles]
        payload: dict[str, Any] = {
            "sub": subject,
            "tenant_id": tenant_id,
            "roles": role_values,
            "email": email,
            "iss": self._settings.jwt_issuer,
            "aud": self._settings.jwt_audience,
            "iat": int(now.timestamp()),
            "exp": int((now + timedelta(minutes=expires_minutes)).timestamp()),
        }
        token = jwt.encode(
            payload,
            self._settings.jwt_secret,
            algorithm=self._settings.jwt_algorithm,
        )
        return str(token)

    def verify(self, token: str) -> Principal:
        try:
            payload = jwt.decode(
                token,
                self._settings.jwt_secret,
                algorithms=[self._settings.jwt_algorithm],
                audience=self._settings.jwt_audience,
                issuer=self._settings.jwt_issuer,
            )
        except JWTError as exc:
            raise AuthenticationError("Invalid or expired token") from exc

        roles_raw = payload.get("roles") or []
        try:
            roles = [Role(r) for r in roles_raw]
        except ValueError as exc:
            raise AuthenticationError("Token contains unknown roles") from exc

        tenant_id = payload.get("tenant_id")
        subject = payload.get("sub")
        if not tenant_id or not subject:
            raise AuthenticationError("Token missing required claims")

        return Principal(
            subject=str(subject),
            tenant_id=str(tenant_id),
            roles=roles,
            email=payload.get("email"),
        )
