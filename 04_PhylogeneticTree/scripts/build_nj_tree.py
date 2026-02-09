# PhylogeneticTree/scripts/build_nj_tree.py
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt
import seaborn as sns

def build_neighbour_joining_trees():
    """Build Neighbour-Joining trees using different distance matrices"""
    
    print("Building Neighbour-Joining trees...")
    
    # Read distance matrices
    euclidean_df = pd.read_csv('output/euclidean_distance_matrix.csv', index_col=0)
    hamming_df = pd.read_csv('output/hamming_distance_matrix.csv', index_col=0)
    
    sample_ids = euclidean_df.index.tolist()
    
    # Convert to BioPython DistanceMatrix format
    def df_to_distance_matrix(df):
        matrix = []
        for i in range(len(df)):
            matrix.append(df.iloc[i, :i+1].tolist())
        return DistanceMatrix(names=df.index.tolist(), matrix=matrix)
    
    # Build trees
    constructor = DistanceTreeConstructor()
    
    # Euclidean distance NJ tree
    print("Building Euclidean distance NJ tree...")
    euclidean_dm = df_to_distance_matrix(euclidean_df)
    nj_tree_euclidean = constructor.nj(euclidean_dm)
    
    # Hamming distance NJ tree
    print("Building Hamming distance NJ tree...")
    hamming_dm = df_to_distance_matrix(hamming_df)
    nj_tree_hamming = constructor.nj(hamming_dm)
    
    # Save trees
    Phylo.write(nj_tree_euclidean, 'output/nj_tree_euclidean.newick', 'newick')
    Phylo.write(nj_tree_hamming, 'output/nj_tree_hamming.newick', 'newick')
    
    # Create visualizations
    plot_nj_trees(nj_tree_euclidean, nj_tree_hamming, sample_ids)
    
    print("Neighbour-Joining trees built and saved")
    return nj_tree_euclidean, nj_tree_hamming

def plot_nj_trees(tree_euclidean, tree_hamming, sample_ids):
    """Plot both NJ trees side by side"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10))
    
    # Plot Euclidean NJ tree
    Phylo.draw(tree_euclidean, axes=ax1, do_show=False)
    ax1.set_title('Neighbour-Joining Tree\n(Euclidean Distance)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    # Plot Hamming NJ tree
    Phylo.draw(tree_hamming, axes=ax2, do_show=False)
    ax2.set_title('Neighbour-Joining Tree\n(Hamming Distance)', 
                 fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig('plots/nj_trees_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('plots/nj_trees_comparison.pdf', bbox_inches='tight')
    
    # Create individual high-quality plots
    fig, ax = plt.subplots(figsize=(15, 10))
    Phylo.draw(tree_euclidean, axes=ax, do_show=False)
    ax.set_title('Neighbour-Joining Tree - Euclidean Distance', 
                fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('plots/nj_tree_euclidean.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    fig, ax = plt.subplots(figsize=(15, 10))
    Phylo.draw(tree_hamming, axes=ax, do_show=False)
    ax.set_title('Neighbour-Joining Tree - Hamming Distance', 
                fontsize=16, fontweight='bold', pad=20)
    plt.tight_layout()
    plt.savefig('plots/nj_tree_hamming.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print("NJ tree plots saved")

if __name__ == "__main__":
    build_neighbour_joining_trees()