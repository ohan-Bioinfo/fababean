#!/bin/bash
# PhylogeneticTree/scripts/run_phylogenetic_analysis.sh

echo "=== Faba Bean Phylogenetic Analysis Pipeline ==="

# Step 1: Convert PLINK to VCF for phylogenetic analysis
echo "Step 1: Converting PLINK to VCF format..."
plink --bfile data/faba_fingerprint \
      --recode vcf \
      --out data/faba_fingerprint --allow-extra-chr

# Step 2: Convert VCF to PHYLIP format
echo "Step 2: Converting to PHYLIP format..."
python scripts/convert_to_phylip.py

# Step 3: Generate distance matrix
echo "Step 3: Generating distance matrix..."
python scripts/generate_distance_matrix.py

# Step 4: Build Neighbour-Joining tree
echo "Step 4: Building Neighbour-Joining tree..."
python scripts/build_nj_tree.py

# Step 5: Build Maximum Likelihood tree
echo "Step 5: Building Maximum Likelihood tree..."
python scripts/build_ml_tree.py

# Step 6: Root trees and create final visualizations
echo "Step 6: Creating rooted trees and visualizations..."
python scripts/root_and_visualize_trees.py

echo "=== Phylogenetic Analysis Complete ==="
echo "Check output/ and plots/ directories for results"