#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

require_cmd plink
require_cmd admixture
require_cmd python3

cd "$PROJECT_ROOT"

base_bfile="01_Raw/03_LD_Prune/Faba_chrOnly_pruned"
if [[ ! -f "${base_bfile}.bed" || ! -f "${base_bfile}.bim" || ! -f "${base_bfile}.fam" ]]; then
  log "ERROR: required LD-pruned files not found under ${base_bfile}.*"
  exit 1
fi

run mkdir -p 02_Diversity/PCA 02_Diversity/IBS 02_Diversity/IBD 02_Diversity/Kinship 02_Diversity/Admixture

log "Stage 02: PCA"
run plink \
  --bfile "$base_bfile" \
  --allow-extra-chr \
  --pca 10 \
  --out 02_Diversity/PCA/Faba_PCA \
  --threads "$THREADS"

log "Stage 02: IBS matrix"
run plink \
  --bfile "$base_bfile" \
  --allow-extra-chr \
  --distance square ibs \
  --out 02_Diversity/IBS/Faba_IBS \
  --threads "$THREADS"

log "Stage 02: IBD/PI-HAT (genome)"
run plink \
  --bfile "$base_bfile" \
  --allow-extra-chr \
  --genome \
  --out 02_Diversity/IBD/Faba_IBD \
  --threads "$THREADS"

if command -v plink2 >/dev/null 2>&1; then
  log "Stage 02: Kinship (KING triangle)"
  run plink2 \
    --bfile "$base_bfile" \
    --make-king triangle \
    --out 02_Diversity/Kinship/Faba_Kinship \
    --threads "$THREADS"
else
  log "WARNING: plink2 not found; skipping KING kinship step."
fi

log "Stage 02: Preparing ADMIXTURE input"
run cp -f "${base_bfile}.bed" 02_Diversity/Admixture/Faba_chrOnly_pruned.bed
run cp -f "${base_bfile}.bim" 02_Diversity/Admixture/Faba_chrOnly_pruned.bim
run cp -f "${base_bfile}.fam" 02_Diversity/Admixture/Faba_chrOnly_pruned.fam
run cp -f 02_Diversity/Admixture/Faba_chrOnly_pruned.bim 02_Diversity/Admixture/Faba_chrOnly_pruned_original.bim

run python3 workflow/scripts/remap_bim_chromosomes.py \
  --mapping 02_Diversity/Admixture/chrom_mapping.txt \
  --input-bim 02_Diversity/Admixture/Faba_chrOnly_pruned_original.bim \
  --output-bim 02_Diversity/Admixture/Faba_chrOnly_pruned_numeric.bim

run cp -f 02_Diversity/Admixture/Faba_chrOnly_pruned.fam 02_Diversity/Admixture/Faba_chrOnly_pruned_numeric.fam
run ln -sfn Faba_chrOnly_pruned.bed 02_Diversity/Admixture/Faba_chrOnly_pruned_numeric.bed

K_MIN="${K_MIN:-2}"
K_MAX="${K_MAX:-10}"
ADMIXTURE_SEED="${ADMIXTURE_SEED:-43}"

log "Stage 02: ADMIXTURE run K=${K_MIN}..${K_MAX}"
if [[ "$DRY_RUN" == "1" ]]; then
  for ((k = K_MIN; k <= K_MAX; k++)); do
    printf '[DRY RUN] (cd 02_Diversity/Admixture && admixture --cv=5 -s %q Faba_chrOnly_pruned_numeric.bed %q > log_K%q.out)\n' "$ADMIXTURE_SEED" "$k" "$k"
  done
else
  for ((k = K_MIN; k <= K_MAX; k++)); do
    log "Running ADMIXTURE K=${k}"
    (
      cd 02_Diversity/Admixture
      admixture --cv=5 -s "$ADMIXTURE_SEED" Faba_chrOnly_pruned_numeric.bed "$k" | tee "log_K${k}.out"
    )
  done
fi

run python3 workflow/scripts/extract_admixture_cv.py \
  --input-dir 02_Diversity/Admixture \
  --output 02_Diversity/Admixture/cv_summary.txt \
  --detailed-output 02_Diversity/Admixture/cv_summary_detailed.tsv

if command -v Rscript >/dev/null 2>&1; then
  if [[ "${RUN_R_PLOTS:-1}" == "1" ]]; then
    log "Stage 02: ADMIXTURE plotting"
    run_in_dir "02_Diversity/Admixture" Rscript visualize_cv.R
  fi
else
  log "WARNING: Rscript not found; skipping ADMIXTURE plotting."
fi

log "Stage 02 complete"
