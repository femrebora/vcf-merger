from phelix_vault.infrastructure.logging.safe import redact_value, safe_extra


def test_redacts_sensitive_keys() -> None:
    assert redact_value("vcf_path", "/data/x.vcf") == "[REDACTED]"
    assert redact_value("password", "secret") == "[REDACTED]"
    extra = safe_extra(asset_id="PHX-AST-1", token="abc")
    assert extra["asset_id"] == "PHX-AST-1"
    assert extra["token"] == "[REDACTED]"
