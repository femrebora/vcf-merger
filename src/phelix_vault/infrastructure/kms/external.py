"""External KMS stub — future AWS KMS / GCP KMS / HashiCorp Vault Transit."""

from __future__ import annotations

from phelix_vault.infrastructure.kms.base import DataKey


class ExternalKMSProvider:
    """Interface placeholder. Do not fake unfinished security features."""

    def __init__(self, *, provider_name: str = "external-kms") -> None:
        self.provider_name = provider_name

    def create_data_key(self, *, key_id: str | None = None) -> DataKey:
        raise NotImplementedError(
            "ExternalKMSProvider is a future integration stub. "
            "Configure kms_provider=development for local use."
        )

    def decrypt_data_key(self, *, key_id: str, encrypted_key: bytes) -> bytes:
        raise NotImplementedError("ExternalKMSProvider is a future integration stub.")

    def rotate_key(self, key_id: str) -> str:
        raise NotImplementedError("ExternalKMSProvider is a future integration stub.")
