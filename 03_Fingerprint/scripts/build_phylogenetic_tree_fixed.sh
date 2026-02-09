#!/bin/bash

# Create necessary directories
mkdir -p output plots

echo "=== Building Phylogenetic Tree using 150 SNPs ==="

# Step 1: Convert VCF to PHYLIP format
echo "Step 1: Converting VCF to PHYLIP format..."
python3 scripts/vcf_to_phylip.py output/faba_150.vcf output/faba_150.phy

# Step 2: Build tree using different methods
echo "Step 2: Building phylogenetic tree..."

# Method 1: Using FastTree (if available)
if command -v FastTree &> /dev/null; then
    echo "Method 1: Building tree with FastTree..."
    FastTree -nt -gtr output/faba_150.phy > output/faba_150_fasttree.nwk
    echo "FastTree analysis completed: output/faba_150_fasttree.nwk"
else
    echo "FastTree not found, skipping FastTree method"
fi

# Method 2: Using R for phylogenetic analysis
echo "Method 2: Building tree with R..."
Rscript scripts/build_tree.R

# Method 3: Using Python for phylogenetic analysis
echo "Method 3: Building tree with Python..."
python3 scripts/simple_tree.py

# Step 4: Compare the trees
echo "Step 4: Comparing generated trees..."
echo "Generated tree files:"
ls -la output/faba_150*.nwk

echo ""
echo "Tree file contents:"
for tree_file in output/faba_150*.nwk; do
    echo "=== $(basename $tree_file) ==="
    cat $tree_file
    echo ""
done

# Step 5: Generate summary
echo "=== Phylogenetic Tree Construction Completed ==="
echo "Output files:"
echo "  - PHYLIP format: output/faba_150.phy"
echo "  - FASTA format: output/faba_150.fasta"
echo "  - FastTree: output/faba_150_fasttree.nwk"
echo "  - Neighbor-Joining tree: output/faba_150_nj_tree.nwk"
echo "  - Hierarchical clustering tree: output/faba_150_hclust_tree.nwk"
echo "  - Python tree: output/faba_150_python_tree.nwk"
echo "  - Plots: plots/faba_150_phylogenetic_tree.pdf"
echo "  - Plots: plots/faba_150_hclust_tree.pdf"
echo "  - Plots: plots/faba_150_python_tree.pdf"

