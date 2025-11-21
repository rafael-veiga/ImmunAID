#!/usr/bin/env bash
set -euo pipefail

# --- helpers ---
log() { echo -e "\033[1;34m[INFO]\033[0m $*"; }
err() { echo -e "\033[1;31m[ERROR]\033[0m $*" >&2; }

# --- paths ---
SCRIPT_DIR="$(cd -- "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd)"
CONTAINERS_DIR="${SCRIPT_DIR}/containers"

# --- docker checks ---
if ! command -v docker >/dev/null 2>&1; then
  err "Docker is not installed."
  exit 1
fi
if ! docker info >/dev/null 2>&1; then
  err "Docker daemon is not running or not accessible. Try: sudo systemctl start docker"
  exit 1
fi

# --- build images ---
log "Building r-container (target: base)..."
docker build -t r-container --target base "${CONTAINERS_DIR}/r_container"

log "Building py-container..."
docker build -t py-container "${CONTAINERS_DIR}/python_container"

log "All images built successfully! 🎉"
