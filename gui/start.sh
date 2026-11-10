#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
if [[ "$(uname -s)" != Linux || "$(uname -m)" != x86_64 ]]; then
  printf '%s\n' 'This first guided launcher supports Linux x86_64 and Windows through WSL2.' 'See gui/README.md for manual installation and platform limitations.'
  exit 1
fi
trap 'result=$?; if [[ $result != 0 && -t 0 ]]; then read -r -p "Setup stopped. Press Enter to close this window." _; fi' EXIT
CONDA_BIN="${CONDA_EXE:-}"
if [[ ! -x "$CONDA_BIN" ]]; then CONDA_BIN="$(command -v conda || true)"; fi
if [[ ! -x "$CONDA_BIN" ]]; then
  for candidate in "$HOME/miniforge3/bin/conda" "$HOME/miniconda3/bin/conda" "$HOME/anaconda3/bin/conda" "$ROOT/.cache/bootstrap/miniforge/bin/conda"; do
    if [[ -x "$candidate" ]]; then CONDA_BIN="$candidate"; break; fi
  done
fi
if [[ ! -x "$CONDA_BIN" ]]; then
  printf '%s\n' 'Conda was not found. The guided setup can install Miniforge inside this project.' 'The first setup downloads the prediction environment and needs several GB of free space.'
  read -r -p 'Install the required environment now? [y/N] ' answer
  if [[ ! "$answer" =~ ^[yY]$ ]]; then exit 1; fi
  BOOTSTRAP="$ROOT/.cache/bootstrap";mkdir -p "$BOOTSTRAP"
  MINIFORGE_VERSION='26.9.4-2'
  INSTALLER="Miniforge3-${MINIFORGE_VERSION}-Linux-x86_64.sh"
  fetch() {
    if command -v curl >/dev/null 2>&1; then curl --fail --location --retry 3 "$1" --output "$2"
    elif command -v wget >/dev/null 2>&1; then wget --tries=3 "$1" -O "$2"
    else printf '%s\n' 'Install curl or wget, then run this launcher again.';return 1;fi
  }
  BASE="https://github.com/conda-forge/miniforge/releases/download/${MINIFORGE_VERSION}"
  fetch "$BASE/$INSTALLER" "$BOOTSTRAP/$INSTALLER"
  fetch "$BASE/$INSTALLER.sha256" "$BOOTSTRAP/$INSTALLER.sha256"
  HASH="$(awk '{print $1;exit}' "$BOOTSTRAP/$INSTALLER.sha256")"
  [[ "$HASH" =~ ^[a-fA-F0-9]{64}$ ]] || { printf '%s\n' 'Invalid installer checksum.';exit 1; }
  printf '%s  %s\n' "$HASH" "$BOOTSTRAP/$INSTALLER" | sha256sum --check
  bash "$BOOTSTRAP/$INSTALLER" -b -p "$BOOTSTRAP/miniforge"
  CONDA_BIN="$BOOTSTRAP/miniforge/bin/conda"
fi
PREFIX="${EVOPHI_APP_PREFIX:-$ROOT/.cache/app-env}"
SPEC_HASH="$(sha256sum "$ROOT/gui/environment.yml" | awk '{print $1}')"
CURRENT_HASH="$(cat "$PREFIX/.evophi-spec" 2>/dev/null || true)"
if [[ ! -x "$PREFIX/bin/python" ]]; then
  printf '%s\n' 'Installing the EvoPhi-DDG prediction environment. This is needed only on first launch.'
  "$CONDA_BIN" env create --prefix "$PREFIX" --file "$ROOT/gui/environment.yml" --yes
  printf '%s\n' "$SPEC_HASH" > "$PREFIX/.evophi-spec"
elif [[ "$CURRENT_HASH" != "$SPEC_HASH" ]]; then
  printf '%s\n' 'Updating the app environment for this version.'
  "$CONDA_BIN" env update --prefix "$PREFIX" --file "$ROOT/gui/environment.yml" --prune
  printf '%s\n' "$SPEC_HASH" > "$PREFIX/.evophi-spec"
fi
export PATH="$PREFIX/bin:$PATH"
exec "$PREFIX/bin/python" "$ROOT/gui/launch.py" "$@"
