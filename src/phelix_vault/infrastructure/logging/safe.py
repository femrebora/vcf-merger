"""Safe logging utilities — never log genomic contents or identity payloads."""

from __future__ import annotations

import logging
import re
from typing import Any
from uuid import uuid4

_SENSITIVE_KEYS = re.compile(
    r"(password|secret|token|authorization|fastq|bam|cram|vcf|sequence|tc_no|national_id)",
    re.IGNORECASE,
)


def new_request_id() -> str:
    return uuid4().hex


def redact_value(key: str, value: Any) -> Any:
    if _SENSITIVE_KEYS.search(key):
        return "[REDACTED]"
    if isinstance(value, str) and len(value) > 500:
        return value[:64] + "...[truncated]"
    return value


def safe_extra(**kwargs: Any) -> dict[str, Any]:
    return {k: redact_value(k, v) for k, v in kwargs.items()}


def get_logger(name: str) -> logging.Logger:
    return logging.getLogger(name)


def configure_logging(level: str = "INFO") -> None:
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s request_id=%(request_id)s %(message)s",
    )
    logging.getLogger().addFilter(_RequestIdFilter())


class _RequestIdFilter(logging.Filter):
    def filter(self, record: logging.LogRecord) -> bool:
        if not hasattr(record, "request_id"):
            record.request_id = "-"
        return True
