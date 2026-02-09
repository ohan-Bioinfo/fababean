# scripts/create_ultraclean_trees.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform
import os

def create_ultraclean_trees():
    """Create ultra-clean, publication-ready trees"""
    
    print("=== CREATING ULTRA-CLEAN PUBLICATION TREES ===")
    
    # Load data
    clusters_df = pd.read_csv('output/phylogenetic_clusters.csv')
    samples = clusters_df['Sample'].tolist()
    distance_df = pd.read_csv('output/genetic_distance_matrix.csv', index_col=0)
    distance_matrix = distance_df.loc[samples, samples].values
    
    # Create linkage
    condensed_dist = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed_dist, method='average')
    
    # Style 1: Nature-style minimal
    plt.figure(figsize=(8, 6))
    dendrogram(linkage_matrix, 
               labels=samples,
               orientation='right',
               leaf_font_size=8,
               leaf_rotation=0,
               color_threshold=0,
               above_threshold_color='black')
    
    plt.xlabel('Genetic Distance', fontsize=9)
    plt.gca().spines['top'].set_visible(False)
    plt.gca().spines['right'].set_visible(False)
    plt.gca().spines['left'].set_visible(False)
    plt.tight_layout()
    plt.savefig('plots/tree_nature_style.png', dpi=400, bbox_inches='tight')
    plt.savefig('plots/tree_nature_style.pdf', dpi=400, bbox_inches='tight')
    plt.show()
    
    # Style 2: High-contrast for presentations
    plt.figure(figsize=(10, 6))
    dendrogram(linkage_matrix, 
               labels=samples,
               orientation='top',
               leaf_rotation=45,
               leaf_font_size=10,
               color_threshold=0.7 * max(linkage_matrix[:, 2]))
    
    plt.ylabel('Genetic Distance', fontsize=11, fontweight='bold')
    plt.title('', fontsize=12)  # No title for clean look
    plt.grid(True, alpha=0.2, axis='y')
    plt.tight_layout()
    plt.savefig('plots/tree_high_contrast.png', dpi=400, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.show()

if __name__ == "__main__":
    create_ultraclean_trees()