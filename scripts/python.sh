#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.."
if [[ -f .runtime/env.sh ]]; then
    source .runtime/env.sh
fi
PYTHON="$(readlink -f .runtime/venv)/bin/python"
exec "$PYTHON" "$@"
