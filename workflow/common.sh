#!/usr/bin/env bash
set -euo pipefail

PROJECT_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THREADS="${THREADS:-8}"
DRY_RUN="${DRY_RUN:-0}"

log() {
  printf '[%s] %s\n' "$(date '+%F %T')" "$*"
}

require_cmd() {
  local cmd="$1"
  if ! command -v "$cmd" >/dev/null 2>&1; then
    if [[ "$DRY_RUN" == "1" ]]; then
      log "WARNING: command not found in DRY_RUN mode: $cmd"
      return 0
    fi
    log "ERROR: required command not found: $cmd"
    exit 1
  fi
}

run() {
  if [[ "$DRY_RUN" == "1" ]]; then
    printf '[DRY RUN] '
    printf '%q ' "$@"
    printf '\n'
  else
    "$@"
  fi
}

run_in_dir() {
  local rel_dir="$1"
  shift
  if [[ "$DRY_RUN" == "1" ]]; then
    printf '[DRY RUN] (cd %s && ' "$rel_dir"
    printf '%q ' "$@"
    printf ')\n'
  else
    (
      cd "$PROJECT_ROOT/$rel_dir"
      "$@"
    )
  fi
}
