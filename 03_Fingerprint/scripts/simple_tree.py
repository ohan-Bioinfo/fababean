#!/usr/bin/env python3
"""
Simple phylogenetic tree construction from genotype data
"""

import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import pdist, squareform
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def main():
    print("Reading genotype data...")
    df = pd.read_csv("output/faba_150_gt.tsv", sep="\t", index_col=0)
    
    print(f"Data shape: {df.shape}")
    print(f"Samples: {list(df.index)}")
    
    # Convert genotypes to numeric
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
    numeric_df = df.map(genotype_to_numeric)
    
    # Fill missing values with mode for each SNP
    numeric_df = numeric_df.apply(lambda x: x.fillna(x.mode()[0] if len(x.mode()) > 0 else 0))
    
    # Calculate distance matrix
    print("Calculating distance matrix...")
    dist_matrix = pdist(numeric_df, metric='hamming')
    dist_square = squareform(dist_matrix)
    
    # Build tree using hierarchical clustering
    print("Building tree...")
    Z = linkage(dist_matrix, method='average')
    
    # Convert to Newick format
    def to_newick(node, parent_dist, leaf_names, is_leaf=False):
        if is_leaf:
            return leaf_names[node]
        else:
            left = to_newick(Z[node, 0], Z[node, 2], leaf_names, Z[node, 0] < len(leaf_names))
            right = to_newick(Z[node, 1], Z[node, 2], leaf_names, Z[node, 1] < len(leaf_names))
            return f"({left}:{Z[node, 2]-parent_dist:.6f},{right}:{Z[node, 2]-parent_dist:.6f})"
    
    tree_newick = to_newick(Z[-1, 0], Z[-1, 2], list(numeric_df.index), Z[-1, 0] < len(numeric_df.index)) + ";"
    
    # Save tree
    with open("output/faba_150_python_tree.nwk", "w") as f:
        f.write(tree_newick)
    
    # Plot dendrogram
    plt.figure(figsize=(12, 8))
    dendrogram(Z, labels=list(numeric_df.index), orientation='right')
    plt.title("Phylogenetic Tree - 150 SNPs (Python)")
    plt.tight_layout()
    plt.savefig("plots/faba_150_python_tree.pdf")
    plt.close()
    
    print("Python tree saved to: output/faba_150_python_tree.nwk")
    print("Python tree plot saved to: plots/faba_150_python_tree.pdf")

if __name__ == "__main__":
    main()
