#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

require_cmd plink
require_cmd python3

cd "$PROJECT_ROOT"

if [[ ! -f 01_Raw/03_LD_Prune/Faba_chrOnly_pruned.bed ]]; then
  log "ERROR: missing input from Stage 01: 01_Raw/03_LD_Prune/Faba_chrOnly_pruned.bed"
  exit 1
fi

if [[ ! -f 03_Fingerprint/scripts/megaScriptall.sh ]]; then
  log "ERROR: missing fingerprint runner script."
  exit 1
fi

run mkdir -p 03_Fingerprint/data 03_Fingerprint/output 03_Fingerprint/plots

log "Stage 03: Fingerprint pipeline"
run_in_dir "03_Fingerprint" bash scripts/megaScriptall.sh

log "Stage 03 complete"
