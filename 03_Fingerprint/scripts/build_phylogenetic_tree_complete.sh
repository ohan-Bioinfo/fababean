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
/usr/bin/Rscript << 'EOF'
# Load required packages
if (!require("ape", quietly = TRUE)) {
    install.packages("ape", repos="http://cran.r-project.org")
}
if (!require("phangorn", quietly = TRUE)) {
    install.packages("phangorn", repos="http://cran.r-project.org")
}
if (!require("adegenet", quietly = TRUE)) {
    install.packages("adegenet", repos="http://cran.r-project.org")
}

library(ape)
library(phangorn)
library(adegenet)

cat("Reading genotype data...\n")

# Read the genotype TSV file
genotype_data <- read.table("output/faba_150_gt.tsv", header=TRUE, row.names=1, sep="\t", na.strings="NA")

# Convert genotypes to numeric (simple approach)
# A/A -> 0, A/T -> 1, T/T -> 2, etc.
convert_genotypes <- function(gt_matrix) {
    num_matrix <- matrix(NA, nrow=nrow(gt_matrix), ncol=ncol(gt_matrix))
    
    for(i in 1:nrow(gt_matrix)) {
        for(j in 1:ncol(gt_matrix)) {
            gt <- as.character(gt_matrix[i,j])
            if(is.na(gt) || gt == "NA") {
                num_matrix[i,j] <- NA
            } else if(gt == "A/A" || gt == "T/T" || gt == "C/C" || gt == "G/G") {
                num_matrix[i,j] <- 0  # Homozygous
            } else if(gt == "A/T" || gt == "T/A" || gt == "A/C" || gt == "C/A" || 
                      gt == "A/G" || gt == "G/A" || gt == "T/C" || gt == "C/T" ||
                      gt == "T/G" || gt == "G/T" || gt == "C/G" || gt == "G/C") {
                num_matrix[i,j] <- 1  # Heterozygous
            } else {
                num_matrix[i,j] <- NA
            }
        }
    }
    return(num_matrix)
}

cat("Converting genotypes to numeric format...\n")
numeric_genotypes <- convert_genotypes(genotype_data)
rownames(numeric_genotypes) <- rownames(genotype_data)

# Calculate distance matrix
cat("Calculating distance matrix...\n")
dist_matrix <- dist(numeric_genotypes, method="manhattan")

# Build Neighbor-Joining tree
cat("Building Neighbor-Joining tree...\n")
nj_tree <- nj(dist_matrix)

# Root the tree (using first sample as outgroup)
rooted_tree <- root(nj_tree, out=1)

# Save the tree
write.tree(rooted_tree, "output/faba_150_nj_tree.nwk")

# Plot the tree
pdf("plots/faba_150_phylogenetic_tree.pdf", width=10, height=8)
plot(rooted_tree, main="Phylogenetic Tree - 150 SNPs\n(Neighbor-Joining)")
add.scale.bar()
dev.off()

# Also try with hclust
cat("Building tree with hierarchical clustering...\n")
hc <- hclust(dist_matrix, method="average")
hc_tree <- as.phylo(hc)

write.tree(hc_tree, "output/faba_150_hclust_tree.nwk")

pdf("plots/faba_150_hclust_tree.pdf", width=10, height=8)
plot(hc_tree, main="Phylogenetic Tree - 150 SNPs\n(Hierarchical Clustering)")
add.scale.bar()
dev.off()

cat("R analysis completed!\n")
cat("Trees saved to: output/faba_150_nj_tree.nwk and output/faba_150_hclust_tree.nwk\n")
EOF

# Step 3: Create a simple Python-based tree as backup
echo "Step 3: Creating backup tree with Python..."
python3 << 'EOF'
import pandas as pd
import numpy as np
from skbio import DistanceMatrix
from skbio.tree import nj
import os

# Read the genotype data
print("Reading genotype data...")
df = pd.read_csv("output/faba_150_gt.tsv", sep="\t", index_col=0)

# Convert genotypes to numeric (simple encoding)
def genotype_to_numeric(gt):
    if pd.isna(gt) or gt == "NA":
        return np.nan
    elif gt in ["A/A", "T/T", "C/C", "G/G"]:
        return 0  # homozygous
    elif gt in ["A/T", "T/A", "A/C", "C/A", "A/G", "G/A", 
                "T/C", "C/T", "T/G", "G/T", "C/G", "G/C"]:
        return 1  # heterozygous
    else:
        return np.nan

print("Converting genotypes...")
numeric_df = df.applymap(genotype_to_numeric)

# Calculate distance matrix
print("Calculating distance matrix...")
from scipy.spatial.distance import pdist, squareform
dist_matrix = pdist(numeric_df.fillna(0), metric='hamming')
dist_square = squareform(dist_matrix)

# Create and save tree using neighbor joining
print("Building tree...")
dm = DistanceMatrix(dist_square, ids=df.index.tolist())
tree = nj(dm)

# Save tree
with open("output/faba_150_python_tree.nwk", "w") as f:
    f.write(str(tree))

print("Python tree saved to: output/faba_150_python_tree.nwk")
EOF

# Step 4: Generate summary
echo "=== Phylogenetic Tree Construction Completed ==="
echo "Output files:"
echo "  - PHYLIP format: output/faba_150.phy"
echo "  - FASTA format: output/faba_150.fasta"
echo "  - Neighbor-Joining tree: output/faba_150_nj_tree.nwk"
echo "  - Hierarchical clustering tree: output/faba_150_hclust_tree.nwk"
echo "  - Python tree: output/faba_150_python_tree.nwk"
echo "  - Plots: plots/faba_150_phylogenetic_tree.pdf"
echo "  - Plots: plots/faba_150_hclust_tree.pdf"

# Display the first tree
if [ -f "output/faba_150_nj_tree.nwk" ]; then
    echo ""
    echo "Newick format tree (first 200 chars):"
    head -c 200 output/faba_150_nj_tree.nwk
    echo ""
fi
