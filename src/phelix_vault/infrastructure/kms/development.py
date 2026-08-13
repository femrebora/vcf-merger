"""Development KMS — deterministic wrapping for local tests only."""

from __future__ import annotations

import hashlib
import hmac
import os

from phelix_vault.infrastructure.kms.base import DataKey


class DevelopmentKeyProvider:
    """NOT for production. Never commit real keys; secret comes from env."""

    def __init__(self, kek_material: str) -> None:
        self._kek = hashlib.sha256(kek_material.encode("utf-8")).digest()
        self._current_key_id = "dev-kek-v1"

    def create_data_key(self, *, key_id: str | None = None) -> DataKey:
        kid = key_id or self._current_key_id
        plaintext = os.urandom(32)
        encrypted = hmac.new(self._kek, plaintext, hashlib.sha256).digest() + plaintext
        return DataKey(key_id=kid, plaintext_key=plaintext, encrypted_key=encrypted)

    def decrypt_data_key(self, *, key_id: str, encrypted_key: bytes) -> bytes:
        del key_id
        if len(encrypted_key) < 32:
            raise ValueError("Invalid encrypted key")
        mac, plaintext = encrypted_key[:32], encrypted_key[32:]
        expected = hmac.new(self._kek, plaintext, hashlib.sha256).digest()
        if not hmac.compare_digest(mac, expected):
            raise ValueError("Invalid encrypted key")
        return plaintext

    def rotate_key(self, key_id: str) -> str:
        # Metadata-ready rotation: new id without rewriting historical ciphertext metadata.
        version = 2
        if key_id.endswith("-v1"):
            return key_id[:-1] + str(version)
        return f"{key_id}-rotated"
