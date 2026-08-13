"""Object storage abstraction."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import PurePosixPath
from typing import BinaryIO, Protocol, runtime_checkable

from phelix_vault.domain.common import DataResidency
from phelix_vault.domain.errors import ValidationError


@dataclass(slots=True, frozen=True)
class StorageLocation:
    provider: str
    bucket: str
    object_key: str
    residency: DataResidency


@runtime_checkable
class ObjectStorage(Protocol):
    provider_name: str
    residency: DataResidency

    def put(
        self, bucket: str, key: str, data: BinaryIO, *, content_type: str | None = None
    ) -> int: ...

    def get(self, bucket: str, key: str) -> BinaryIO: ...

    def delete(self, bucket: str, key: str) -> None: ...

    def exists(self, bucket: str, key: str) -> bool: ...

    def checksum_sha256(self, bucket: str, key: str) -> str: ...

    def signed_download(self, bucket: str, key: str, *, expires_seconds: int) -> str: ...

    def signed_upload(self, bucket: str, key: str, *, expires_seconds: int) -> str: ...


_FORBIDDEN_KEY_PARTS = ("..", "\\", "\x00")


def validate_object_key(key: str) -> str:
    if not key or key.startswith("/") or any(p in key for p in _FORBIDDEN_KEY_PARTS):
        raise ValidationError("Unsafe object key rejected")
    normalized = PurePosixPath(key).as_posix()
    if normalized.startswith("..") or "/../" in f"/{normalized}/":
        raise ValidationError("Unsafe object key rejected")
    return normalized
