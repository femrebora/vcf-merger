#!/usr/bin/env bash
set -euo pipefail

BASE="${BASE_URL:-http://127.0.0.1:8000}"

echo "Health: $BASE/health"
curl -fsS "$BASE/health" | tee /tmp/phelix-health.json
echo
