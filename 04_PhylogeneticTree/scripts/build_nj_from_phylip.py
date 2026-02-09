# scripts/build_nj_from_phylip.py
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt
import os

def build_nj_trees():
    """Build Neighbour-Joining trees from distance matrix"""
    
    print("Building Neighbour-Joining trees...")
    
    # Ensure directories exist
    os.makedirs('output', exist_ok=True)
    os.makedirs('plots', exist_ok=True)
    
    # Read distance matrix
    try:
        dist_df = pd.read_csv('output/hamming_distance_matrix.csv', index_col=0)
        sample_ids = dist_df.index.tolist()
        print(f"✓ Loaded distance matrix for {len(sample_ids)} samples")
    except Exception as e:
        print(f"✗ Error loading distance matrix: {e}")
        return False
    
    # Convert to BioPython DistanceMatrix format
    def df_to_distance_matrix(df):
        matrix = []
        for i in range(len(df)):
            matrix.append(df.iloc[i, :i+1].tolist())
        return DistanceMatrix(names=df.index.tolist(), matrix=matrix)
    
    try:
        # Build NJ tree
        constructor = DistanceTreeConstructor()
        distance_matrix = df_to_distance_matrix(dist_df)
        nj_tree = constructor.nj(distance_matrix)
        
        # Save unrooted tree
        Phylo.write(nj_tree, 'output/nj_tree_unrooted.newick', 'newick')
        print("✓ NJ tree built and saved")
        
        # Create visualization
        plot_tree(nj_tree, 'Neighbour-Joining Tree (Unrooted)', 'plots/nj_tree_unrooted.png')
        
        return nj_tree, sample_ids
        
    except Exception as e:
        print(f"✗ Error building NJ tree: {e}")
        return False

def plot_tree(tree, title, filename):
    """Plot phylogenetic tree"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    
    fig, ax = plt.subplots(figsize=(12, 8))
    
    # Draw the tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    # Customize the plot
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"✓ Tree plot saved: {filename}")

def main():
    result = build_nj_trees()
    if result:
        print("\n✓ Neighbour-Joining analysis completed successfully")
    else:
        print("\n✗ Neighbour-Joining analysis failed")

if __name__ == "__main__":
    main()