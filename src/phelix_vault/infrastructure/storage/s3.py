"""S3-compatible object storage (MinIO / AWS / private S3)."""

from __future__ import annotations

import hashlib
from typing import BinaryIO, cast

import boto3
from botocore.client import Config
from botocore.exceptions import ClientError

from phelix_vault.domain.common import DataResidency
from phelix_vault.domain.errors import NotFoundError
from phelix_vault.infrastructure.storage.base import validate_object_key


class S3ObjectStorage:
    provider_name = "s3"

    def __init__(
        self,
        *,
        bucket: str,
        region: str,
        residency: DataResidency,
        endpoint_url: str | None = None,
        access_key: str | None = None,
        secret_key: str | None = None,
    ) -> None:
        self.default_bucket = bucket
        self.residency = residency
        session_kwargs: dict[str, str] = {}
        if access_key and secret_key:
            session_kwargs["aws_access_key_id"] = access_key
            session_kwargs["aws_secret_access_key"] = secret_key
        self._client = boto3.client(
            "s3",
            region_name=region,
            endpoint_url=endpoint_url,
            config=Config(signature_version="s3v4"),
            **session_kwargs,
        )  # type: ignore[call-overload]

    def put(self, bucket: str, key: str, data: BinaryIO, *, content_type: str | None = None) -> int:
        safe = validate_object_key(key)
        extra = {"ContentType": content_type} if content_type else {}
        # Stream via upload_fileobj — do not load entire object into memory.
        self._client.upload_fileobj(data, bucket, safe, ExtraArgs=extra or None)
        head = self._client.head_object(Bucket=bucket, Key=safe)
        return int(head.get("ContentLength", 0))

    def get(self, bucket: str, key: str) -> BinaryIO:
        safe = validate_object_key(key)
        try:
            obj = self._client.get_object(Bucket=bucket, Key=safe)
        except ClientError as exc:
            raise NotFoundError("Object not found") from exc
        return cast(BinaryIO, obj["Body"])

    def delete(self, bucket: str, key: str) -> None:
        safe = validate_object_key(key)
        self._client.delete_object(Bucket=bucket, Key=safe)

    def exists(self, bucket: str, key: str) -> bool:
        safe = validate_object_key(key)
        try:
            self._client.head_object(Bucket=bucket, Key=safe)
            return True
        except ClientError:
            return False

    def checksum_sha256(self, bucket: str, key: str) -> str:
        digest = hashlib.sha256()
        body = self.get(bucket, key)
        try:
            while True:
                chunk = body.read(1024 * 1024)
                if not chunk:
                    break
                digest.update(chunk)
        finally:
            close = getattr(body, "close", None)
            if callable(close):
                close()
        return digest.hexdigest()

    def signed_download(self, bucket: str, key: str, *, expires_seconds: int) -> str:
        safe = validate_object_key(key)
        url = self._client.generate_presigned_url(
            "get_object",
            Params={"Bucket": bucket, "Key": safe},
            ExpiresIn=expires_seconds,
        )
        return str(url)

    def signed_upload(self, bucket: str, key: str, *, expires_seconds: int) -> str:
        safe = validate_object_key(key)
        url = self._client.generate_presigned_url(
            "put_object",
            Params={"Bucket": bucket, "Key": safe},
            ExpiresIn=expires_seconds,
        )
        return str(url)
