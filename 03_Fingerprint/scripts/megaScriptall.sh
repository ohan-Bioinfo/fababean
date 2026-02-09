#!/usr/bin/env bash
set -euo pipefail

echo "=== Faba Bean Fingerprint Pipeline ==="

# Step 1: QC filtering
echo "Step 1: Applying QC filters..."
plink --bfile ../01_Raw/03_LD_Prune/Faba_chrOnly_pruned \
      --maf 0.3 \
      --geno 0.05 \
      --mind 0.25 \
      --make-bed \
      --out data/Faba_high_quality

# Step 2: Calculate PIC
echo "Step 2: Calculating PIC values..."
python3 scripts/calculate_pic_complete.py

# Step 3: Select top SNPs
echo "Step 3: Selecting top 150 SNPs..."
python3 scripts/select_top_snps.py

# Step 4: Create fingerprint panel
echo "Step 4: Creating fingerprint panel..."
plink --bfile data/Faba_high_quality \
      --extract data/top_150_snps_list.txt \
      --make-bed \
      --out data/faba_fingerprint

# Step 5: Generate genotype tables
echo "Step 5: Generating genotype tables..."
plink --bfile data/faba_fingerprint \
      --recode A \
      --out output/faba_fingerprint_minimal

plink --bfile data/faba_fingerprint \
      --recode \
      --out output/faba_fingerprint_genotypes

# Step 6: Create visualization
echo "Step 6: Creating fingerprint heatmap..."
python3 scripts/create_fingerprint_heatmap.py

echo "=== Pipeline Complete ==="
echo "Fingerprint panel with 150 SNPs created successfully!"
echo "Check output/ and plots/ directories for results."
