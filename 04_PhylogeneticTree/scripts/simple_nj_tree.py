# scripts/simple_nj_tree.py
import pandas as pd
import numpy as np
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio import Phylo
import matplotlib.pyplot as plt
import os

def build_simple_nj_tree():
    """Build a simple NJ tree as backup"""
    
    print("Building simple Neighbour-Joining tree...")
    
    # Ensure directories exist
    os.makedirs('output', exist_ok=True)
    os.makedirs('plots', exist_ok=True)
    
    try:
        # Read distance matrix
        dist_df = pd.read_csv('output/hamming_distance_matrix.csv', index_col=0)
        
        # Use FID only
        sample_ids = [sid.split('_')[0] if '_' in sid else sid for sid in dist_df.index.tolist()]
        
        print(f"Building tree for {len(sample_ids)} samples: {sample_ids}")
        
        # Convert to BioPython DistanceMatrix format
        def df_to_distance_matrix(df):
            matrix = []
            for i in range(len(df)):
                matrix.append(df.iloc[i, :i+1].tolist())
            return DistanceMatrix(names=sample_ids, matrix=matrix)
        
        # Build NJ tree
        constructor = DistanceTreeConstructor()
        distance_matrix = df_to_distance_matrix(dist_df)
        nj_tree = constructor.nj(distance_matrix)
        
        # Save tree
        Phylo.write(nj_tree, 'output/nj_tree_simple.newick', 'newick')
        print("✓ Simple NJ tree saved")
        
        # Create simple visualization
        plot_simple_tree(nj_tree, 'Neighbour-Joining Tree (Simple)', 'plots/nj_tree_simple.png')
        
        # Create rooted version
        rooted_tree = root_tree_simple(nj_tree)
        if rooted_tree:
            Phylo.write(rooted_tree, 'output/nj_tree_rooted_simple.newick', 'newick')
            plot_simple_tree(rooted_tree, 'Rooted Neighbour-Joining Tree', 'plots/nj_tree_rooted_simple.png')
        
        return True
        
    except Exception as e:
        print(f"✗ Error building simple NJ tree: {e}")
        return False

def plot_simple_tree(tree, title, filename):
    """Create a simple tree plot"""
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    
    fig, ax = plt.subplots(figsize=(12, 8))
    Phylo.draw(tree, axes=ax, do_show=False)
    
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    
    # Make labels bold
    for text in ax.texts:
        text.set_fontweight('bold')
        text.set_fontsize(10)
    
    plt.tight_layout()
    plt.savefig(filename, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✓ Simple tree plot saved: {filename}")

def root_tree_simple(tree):
    """Simple midpoint rooting"""
    try:
        terminals = list(tree.get_terminals())
        
        # Find two most distant terminals
        max_dist = 0
        t1, t2 = None, None
        
        for i in range(len(terminals)):
            for j in range(i+1, len(terminals)):
                dist = tree.distance(terminals[i], terminals[j])
                if dist > max_dist:
                    max_dist = dist
                    t1, t2 = terminals[i], terminals[j]
        
        if t1 and t2:
            tree.root_at_midpoint(t1, t2, max_dist/2)
            print("✓ Tree rooted at midpoint")
            return tree
        else:
            print("⚠ Could not root tree - using unrooted")
            return tree
            
    except Exception as e:
        print(f"⚠ Rooting failed: {e} - using unrooted tree")
        return tree

if __name__ == "__main__":
    success = build_simple_nj_tree()
    if success:
        print("\n✓ Simple NJ tree analysis completed")
    else:
        print("\n✗ Simple NJ tree analysis failed")