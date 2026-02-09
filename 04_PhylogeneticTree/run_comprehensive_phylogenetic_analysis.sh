#!/bin/bash
# scripts/run_comprehensive_phylogenetic_analysis.sh

echo "=== Comprehensive Faba Bean Phylogenetic Analysis ==="

# Step 1: Convert VCF to PHYLIP
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
    python scripts/convert_raw_to_phylip.py
    if [ ! -f "data/faba_fingerprint.phy" ]; then
        exit 1
    fi
fi

# Step 2: Generate distance matrix
echo "Step 2: Generating distance matrix..."
python scripts/generate_distance_from_phylip.py

# Step 3: Build Enhanced NJ trees (both unrooted and rooted)
echo "Step 3: Building Enhanced Neighbour-Joining trees..."
python scripts/build_enhanced_nj_trees.py

# Step 4: Build Enhanced ML trees (both unrooted and rooted)
echo "Step 4: Building Enhanced Maximum Likelihood trees..."
python scripts/build_enhanced_ml_trees.py

# Step 5: Create comprehensive comparison
echo "Step 5: Creating comprehensive tree comparison..."
python scripts/create_comprehensive_tree_comparison.py

echo ""
echo "=== COMPREHENSIVE PHYLOGENETIC ANALYSIS COMPLETE ==="
echo "✓ All tree types generated:"
echo "  - Neighbour-Joining (unrooted and rooted)"
echo "  - Maximum Likelihood (multiple models)"
echo "  - Comprehensive comparisons"
echo ""
echo "✓ Output files:"
echo "  - Tree files: output/*.newick"
echo "  - Visualizations: plots/*.png"
echo "  - Summary: output/phylogenetic_trees_summary.csv"
echo ""
echo "✓ Features:"
echo "  - FID-only sample labels"
echo "  - Bold and larger fonts"
echo "  - Enhanced visualizations"
echo "  - Multiple tree layouts"
echo "  - Publication-ready figures"