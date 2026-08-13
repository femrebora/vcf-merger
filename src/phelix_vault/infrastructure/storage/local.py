"""Local filesystem object storage for development/testing."""

from __future__ import annotations

import hashlib
import shutil
from io import BytesIO
from pathlib import Path
from typing import BinaryIO
from urllib.parse import quote

from phelix_vault.domain.common import DataResidency
from phelix_vault.domain.errors import NotFoundError, ValidationError
from phelix_vault.infrastructure.storage.base import validate_object_key


class LocalObjectStorage:
    provider_name = "local"

    def __init__(self, root: Path, *, residency: DataResidency = DataResidency.TR) -> None:
        self.root = root
        self.residency = residency
        self.root.mkdir(parents=True, exist_ok=True)

    def _path(self, bucket: str, key: str) -> Path:
        safe_key = validate_object_key(key)
        if "/" in bucket or ".." in bucket or not bucket:
            raise ValidationError("Invalid bucket name")
        path = (self.root / bucket / safe_key).resolve()
        root_resolved = (self.root / bucket).resolve()
        if not str(path).startswith(str(root_resolved)):
            raise ValidationError("Unsafe object key rejected")
        return path

    def put(self, bucket: str, key: str, data: BinaryIO, *, content_type: str | None = None) -> int:
        path = self._path(bucket, key)
        path.parent.mkdir(parents=True, exist_ok=True)
        size = 0
        with path.open("wb") as fh:
            while True:
                chunk = data.read(1024 * 1024)
                if not chunk:
                    break
                fh.write(chunk)
                size += len(chunk)
        if content_type:
            meta = path.with_suffix(path.suffix + ".content-type")
            meta.write_text(content_type, encoding="utf-8")
        return size

    def get(self, bucket: str, key: str) -> BinaryIO:
        path = self._path(bucket, key)
        if not path.is_file():
            raise NotFoundError("Object not found")
        return path.open("rb")

    def delete(self, bucket: str, key: str) -> None:
        path = self._path(bucket, key)
        if path.is_file():
            path.unlink()
        ct = path.with_suffix(path.suffix + ".content-type")
        if ct.is_file():
            ct.unlink()

    def exists(self, bucket: str, key: str) -> bool:
        return self._path(bucket, key).is_file()

    def checksum_sha256(self, bucket: str, key: str) -> str:
        digest = hashlib.sha256()
        with self.get(bucket, key) as fh:
            while True:
                chunk = fh.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
        return digest.hexdigest()

    def signed_download(self, bucket: str, key: str, *, expires_seconds: int) -> str:
        # Local "signed" URL is an opaque vault-relative token path for development.
        safe = validate_object_key(key)
        return f"local://{bucket}/{quote(safe)}?expires={expires_seconds}"

    def signed_upload(self, bucket: str, key: str, *, expires_seconds: int) -> str:
        safe = validate_object_key(key)
        return f"local-upload://{bucket}/{quote(safe)}?expires={expires_seconds}"

    def put_bytes(
        self,
        bucket: str,
        key: str,
        content: bytes,
        *,
        content_type: str | None = None,
    ) -> int:
        return self.put(bucket, key, BytesIO(content), content_type=content_type)

    def wipe_bucket(self, bucket: str) -> None:
        target = self.root / bucket
        if target.exists():
            shutil.rmtree(target)
