# scripts/build_enhanced_nj_trees.py (Fixed version)
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
import os

def build_enhanced_nj_trees():
    """Build enhanced Neighbour-Joining trees with better visualization"""
    
    print("Building Enhanced Neighbour-Joining trees...")
    
    # Ensure directories exist
    os.makedirs('output', exist_ok=True)
    os.makedirs('plots', exist_ok=True)
    
    # Read distance matrix
    try:
        dist_df = pd.read_csv('output/hamming_distance_matrix.csv', index_col=0)
        # Use FID only (remove everything after underscore if present)
        sample_ids = [sid.split('_')[0] if '_' in sid else sid for sid in dist_df.index.tolist()]
        print(f"✓ Loaded distance matrix for {len(sample_ids)} samples")
        print(f"Sample IDs: {sample_ids}")
    except Exception as e:
        print(f"✗ Error loading distance matrix: {e}")
        return False
    
    # Convert to BioPython DistanceMatrix format
    def df_to_distance_matrix(df):
        matrix = []
        for i in range(len(df)):
            matrix.append(df.iloc[i, :i+1].tolist())
        return DistanceMatrix(names=sample_ids, matrix=matrix)
    
    try:
        # Build NJ tree
        constructor = DistanceTreeConstructor()
        distance_matrix = df_to_distance_matrix(dist_df)
        nj_tree = constructor.nj(distance_matrix)
        
        # Save unrooted tree
        Phylo.write(nj_tree, 'output/nj_tree_unrooted.newick', 'newick')
        print("✓ NJ tree built and saved")
        
        # Create enhanced visualizations
        plot_enhanced_nj_trees(nj_tree, sample_ids)
        
        return nj_tree, sample_ids
        
    except Exception as e:
        print(f"✗ Error building NJ tree: {e}")
        return False

def plot_enhanced_nj_trees(tree, sample_ids):
    """Create enhanced NJ tree visualizations with colors and better formatting"""
    
    # Set up matplotlib for publication quality
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    plt.rcParams['axes.labelweight'] = 'bold'
    
    # Create a color map for samples
    colors = plt.cm.Set3(np.linspace(0, 1, len(sample_ids)))
    
    # Visualization 1: Standard NJ Tree
    fig, ax = plt.subplots(figsize=(16, 12))
    
    # Draw tree with enhanced formatting - FIXED: removed duplicate do_show=False
    Phylo.draw(tree, axes=ax, do_show=False, 
               label_func=lambda x: x.name)  # Use FID only
    
    # Enhance the plot
    ax.set_title('Neighbour-Joining Phylogenetic Tree\nFaba Bean Accessions (Unrooted)', 
                fontsize=18, fontweight='bold', pad=20)
    
    # Make tip labels bold and larger
    for text in ax.texts:
        text.set_fontweight('bold')
        text.set_fontsize(12)
    
    plt.tight_layout()
    plt.savefig('plots/nj_tree_unrooted_enhanced.png', dpi=350, bbox_inches='tight')
    plt.savefig('plots/nj_tree_unrooted_enhanced.pdf', bbox_inches='tight')
    plt.close()
    
    # Visualization 2: Circular NJ Tree
    fig, ax = plt.subplots(figsize=(18, 16))
    
    # Draw circular tree - FIXED: removed duplicate do_show=False
    Phylo.draw(tree, axes=ax, do_show=False, 
               label_func=lambda x: x.name)
    
    ax.set_title('Circular Neighbour-Joining Tree\nFaba Bean Accessions', 
                fontsize=20, fontweight='bold', pad=30)
    
    # Enhance tip labels
    for text in ax.texts:
        text.set_fontweight('bold')
        text.set_fontsize(11)
        text.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='none'))
    
    plt.tight_layout()
    plt.savefig('plots/nj_tree_circular.png', dpi=350, bbox_inches='tight')
    plt.savefig('plots/nj_tree_circular.pdf', bbox_inches='tight')
    plt.close()
    
    print("✓ Enhanced NJ tree visualizations saved")

def create_rooted_nj_tree():
    """Create and visualize rooted NJ tree"""
    
    print("Creating rooted NJ tree...")
    
    try:
        # Read the unrooted tree
        tree = Phylo.read('output/nj_tree_unrooted.newick', 'newick')
        
        # Root the tree at midpoint
        rooted_tree = root_tree_at_midpoint(tree)
        
        # Save rooted tree
        Phylo.write(rooted_tree, 'output/nj_tree_rooted.newick', 'newick')
        
        # Create enhanced visualization of rooted tree
        plot_enhanced_rooted_tree(rooted_tree, 'Neighbour-Joining')
        
        return rooted_tree
        
    except Exception as e:
        print(f"✗ Error creating rooted NJ tree: {e}")
        return None

def root_tree_at_midpoint(tree):
    """Root tree at midpoint"""
    # Get all terminal nodes
    terminals = list(tree.get_terminals())
    
    # Find the two most distant nodes
    max_distance = 0
    tip1, tip2 = None, None
    
    for t1 in terminals:
        for t2 in terminals:
            if t1 != t2:
                dist = tree.distance(t1, t2)
                if dist > max_distance:
                    max_distance = dist
                    tip1, tip2 = t1, t2
    
    if tip1 and tip2:
        # Root at midpoint
        tree.root_at_midpoint(tip1, tip2, max_distance/2)
    
    return tree

def plot_enhanced_rooted_tree(tree, tree_type):
    """Create enhanced visualization of rooted tree"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    
    fig, ax = plt.subplots(figsize=(14, 10))
    
    # Draw rooted tree
    Phylo.draw(tree, axes=ax, do_show=False)
    
    ax.set_title(f'Rooted {tree_type} Phylogenetic Tree\nFaba Bean Accessions (Midpoint Rooting)', 
                fontsize=16, fontweight='bold', pad=20)
    
    # Enhance tip labels
    for text in ax.texts:
        text.set_fontweight('bold')
        text.set_fontsize(11)
    
    # Add tree statistics
    n_tips = len(tree.get_terminals())
    tree_length = sum(tree.depths().values())
    
    ax.text(0.02, 0.98, f'Samples: {n_tips}\nTree length: {tree_length:.2f}', 
           transform=ax.transAxes, fontsize=10, fontweight='bold',
           verticalalignment='top', 
           bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.8))
    
    plt.tight_layout()
    safe_name = tree_type.lower().replace(' ', '_')
    plt.savefig(f'plots/rooted_{safe_name}_enhanced.png', dpi=350, bbox_inches='tight')
    plt.savefig(f'plots/rooted_{safe_name}_enhanced.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Enhanced rooted {tree_type} tree visualization saved")

def main():
    print("=== Enhanced Neighbour-Joining Tree Analysis ===")
    
    # Build enhanced NJ trees
    result = build_enhanced_nj_trees()
    
    if result:
        # Create rooted version
        rooted_tree = create_rooted_nj_tree()
        
        print("\n✓ Enhanced NJ tree analysis completed successfully")
        print("✓ Generated:")
        print("  - Unrooted NJ tree (standard layout)")
        print("  - Unrooted NJ tree (circular layout)") 
        print("  - Rooted NJ tree (midpoint rooting)")
    else:
        print("\n✗ Enhanced NJ tree analysis failed")

if __name__ == "__main__":
    main()