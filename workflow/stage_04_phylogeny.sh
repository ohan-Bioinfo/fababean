#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

require_cmd plink
require_cmd python3

cd "$PROJECT_ROOT"

if [[ ! -f 03_Fingerprint/data/faba_fingerprint.bed ]]; then
  log "ERROR: missing fingerprint panel bed file from Stage 03."
  exit 1
fi

run mkdir -p 04_PhylogeneticTree/data 04_PhylogeneticTree/output 04_PhylogeneticTree/plots

log "Stage 04: Sync fingerprint panel into phylogeny input folder"
run cp -f 03_Fingerprint/data/faba_fingerprint.bed 04_PhylogeneticTree/data/
run cp -f 03_Fingerprint/data/faba_fingerprint.bim 04_PhylogeneticTree/data/
run cp -f 03_Fingerprint/data/faba_fingerprint.fam 04_PhylogeneticTree/data/

log "Stage 04: Generate VCF for phylogenetic conversion"
run plink \
  --bfile 04_PhylogeneticTree/data/faba_fingerprint \
  --recode vcf \
  --allow-extra-chr \
  --out 04_PhylogeneticTree/data/faba_fingerprint

log "Stage 04: Comprehensive phylogenetic analysis"
run_in_dir "04_PhylogeneticTree" bash run_comprehensive_phylogenetic_analysis.sh

log "Stage 04 complete"
