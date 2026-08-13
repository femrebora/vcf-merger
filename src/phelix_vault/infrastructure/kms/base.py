"""Key management abstractions (envelope-encryption ready)."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable


@dataclass(slots=True, frozen=True)
class DataKey:
    key_id: str
    plaintext_key: bytes
    encrypted_key: bytes
    algorithm: str = "AES-256-GCM"


@runtime_checkable
class KeyManagementProvider(Protocol):
    def create_data_key(self, *, key_id: str | None = None) -> DataKey: ...

    def decrypt_data_key(self, *, key_id: str, encrypted_key: bytes) -> bytes: ...

    def rotate_key(self, key_id: str) -> str:
        """Return new key id; callers re-wrap DEKs asynchronously."""
        ...
