#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

bash "$SCRIPT_DIR/stage_01_qc.sh"
bash "$SCRIPT_DIR/stage_02_diversity.sh"
bash "$SCRIPT_DIR/stage_03_fingerprint.sh"
bash "$SCRIPT_DIR/stage_04_phylogeny.sh"
