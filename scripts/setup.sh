#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
RUNTIME="${PROTEINLOC_RUNTIME:-${XDG_DATA_HOME:-$HOME/.local/share}/proteinloc}"
mkdir -p "$ROOT/.runtime/bin" "$RUNTIME"
cd "$ROOT"
export PATH="$ROOT/.runtime/bin:$HOME/.local/bin:$PATH"
if ! command -v uv >/dev/null 2>&1; then
    curl --fail --location --retry 2 https://astral.sh/uv/install.sh -o .runtime/install-uv.sh
    UV_INSTALL_DIR="$ROOT/.runtime/bin" UV_NO_MODIFY_PATH=1 sh .runtime/install-uv.sh
fi
if [[ ! -x "$RUNTIME/venv/bin/python" ]]; then
    export UV_PYTHON_INSTALL_DIR="$RUNTIME/python"
    uv python install 3.12.13
    PYTHON="$(uv python find 3.12.13)"
    uv venv --python "$PYTHON" "$RUNTIME/venv"
fi
if [[ -d "$ROOT/.runtime/wheels" ]] && compgen -G "$ROOT/.runtime/wheels/numpy-*.whl" >/dev/null; then
    uv pip install --python "$RUNTIME/venv/bin/python" --no-index --find-links "$ROOT/.runtime/wheels" -r requirements.txt
else
    uv pip install --python "$RUNTIME/venv/bin/python" -r requirements.txt
fi
if [[ ! -e .runtime/venv ]]; then
    ln -s "$RUNTIME/venv" .runtime/venv
fi
.runtime/venv/bin/python -c 'import tkinter, numpy, openpyxl; print("GUI / Excel / PredGPI dependencies ready")'
.runtime/venv/bin/python scripts/prepare_tk.py
.runtime/venv/bin/python scripts/setup_support.py
echo "运行 bash run.sh 启动 GUI；运行 .runtime/venv/bin/python pipeline.py doctor 查看模型状态。"
