#!/usr/bin/env bash
set -euo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/common.sh"

require_cmd plink
require_cmd python3

cd "$PROJECT_ROOT"

input_vcf="00_GBS/pop.vcf.gz"
if [[ ! -f "$input_vcf" ]]; then
  log "ERROR: input VCF not found: $input_vcf"
  exit 1
fi

run mkdir -p 01_Raw/02.2__SNP_Filter 01_Raw/03_LD_Prune

log "Stage 01: VCF -> PLINK conversion"
run plink \
  --vcf "$input_vcf" \
  --snps-only just-acgt \
  --allow-extra-chr \
  --double-id \
  --set-missing-var-ids '@:#:$1:$2' \
  --make-bed \
  --out 01_Raw/Faba_chrOnly_raw \
  --threads "$THREADS"

log "Stage 01: Missingness and heterozygosity reports"
run plink \
  --bfile 01_Raw/Faba_chrOnly_raw \
  --allow-extra-chr \
  --missing \
  --out 01_Raw/02.2__SNP_Filter/Faba_missingness \
  --threads "$THREADS"

run plink \
  --bfile 01_Raw/Faba_chrOnly_raw \
  --allow-extra-chr \
  --het \
  --out 01_Raw/02.2__SNP_Filter/Faba_het \
  --threads "$THREADS"

log "Stage 01: Sample/SNP filtering"
run plink \
  --bfile 01_Raw/Faba_chrOnly_raw \
  --allow-extra-chr \
  --mind 0.50 \
  --make-bed \
  --out 01_Raw/02.2__SNP_Filter/Faba_chrOnly_mind05 \
  --threads "$THREADS"

run plink \
  --bfile 01_Raw/02.2__SNP_Filter/Faba_chrOnly_mind05 \
  --allow-extra-chr \
  --geno 0.05 \
  --maf 0.05 \
  --make-bed \
  --out 01_Raw/02.2__SNP_Filter/Faba_chrOnly_geno05_maf05 \
  --threads "$THREADS"

log "Stage 01: LD pruning"
run plink \
  --bfile 01_Raw/02.2__SNP_Filter/Faba_chrOnly_geno05_maf05 \
  --allow-extra-chr \
  --indep-pairwise 50 5 0.5 \
  --out 01_Raw/03_LD_Prune/Faba_chrOnly_ld \
  --threads "$THREADS"

run plink \
  --bfile 01_Raw/02.2__SNP_Filter/Faba_chrOnly_geno05_maf05 \
  --allow-extra-chr \
  --extract 01_Raw/03_LD_Prune/Faba_chrOnly_ld.prune.in \
  --make-bed \
  --out 01_Raw/03_LD_Prune/Faba_chrOnly_pruned \
  --threads "$THREADS"

log "Stage 01: Optional QC visualizations"
if [[ -f 01_Raw/03_LD_Prune/vsdistribution.py ]]; then
  run_in_dir "01_Raw/03_LD_Prune" python3 vsdistribution.py
fi

if [[ -f 00_GBS/arc/03.snp.stat.xls && -f 00_GBS/arc/04indel.stat.xls ]]; then
  run_in_dir "00_GBS/arc" python3 ../vs.py
fi

log "Stage 01 complete"
