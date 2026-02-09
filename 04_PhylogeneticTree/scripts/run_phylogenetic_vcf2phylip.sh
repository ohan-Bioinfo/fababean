#!/bin/bash
# scripts/run_phylogenetic_vcf2phylip.sh

echo "=== Faba Bean Phylogenetic Analysis (vcf2phylip method) ==="

# Step 1: Convert VCF to PHYLIP using vcf2phylip.py
echo "Step 1: Converting VCF to PHYLIP format..."
python vcf2phylip.py -i data/faba_fingerprint.vcf -o data/faba_fingerprint

# Check if conversion was successful
if [ -f "faba_fingerprint.min4.phy" ]; then
    mv faba_fingerprint.min4.phy data/faba_fingerprint.phy
    echo "✓ PHYLIP file created: data/faba_fingerprint.phy"
elif [ -f "data/faba_fingerprint.phy" ]; then
    echo "✓ PHYLIP file exists: data/faba_fingerprint.phy"
else
    echo "ERROR: vcf2phylip conversion failed"
    # Try alternative conversion method
    echo "Trying alternative conversion..."
    python scripts/convert_raw_to_phylip.py
    if [ ! -f "data/faba_fingerprint.phy" ]; then
        exit 1
    fi
fi

# Verify PHYLIP file
echo "PHYLIP file contents:"
head -n 3 data/faba_fingerprint.phy

# Step 2: Generate distance matrix from PHYLIP
echo "Step 2: Generating distance matrix..."
python scripts/generate_distance_from_phylip.py

# Step 3: Build Neighbour-Joining tree
echo "Step 3: Building Neighbour-Joining tree..."
python scripts/build_nj_from_phylip.py

# Step 4: Build Maximum Likelihood tree
echo "Step 4: Building Maximum Likelihood tree..."
python scripts/build_ml_from_phylip.py

# Step 5: Root trees and create visualizations
echo "Step 5: Creating rooted trees and visualizations..."
python scripts/root_and_visualize_from_phylip.py

echo "=== Phylogenetic Analysis Complete ==="
echo "Check output/ and plots/ directories for results"