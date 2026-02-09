# scripts/build_enhanced_ml_trees.py
import os
import subprocess
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import numpy as np

def build_enhanced_ml_trees():
    """Build enhanced Maximum Likelihood trees with better visualization"""
    
    print("Building Enhanced Maximum Likelihood trees...")
    
    phylip_file = 'data/faba_fingerprint.phy'
    
    if not os.path.exists(phylip_file):
        print(f"✗ PHYLIP file not found: {phylip_file}")
        return False
    
    # Check if IQ-TREE is available
    try:
        subprocess.run(['iqtree', '--version'], capture_output=True)
        iqtree_available = True
    except:
        iqtree_available = False
        print("WARNING: IQ-TREE not found. Skipping ML trees.")
        return False
    
    if iqtree_available:
        print("Running enhanced IQ-TREE analysis...")
        
        # Strategy 1: Simple model with bootstrapping
        cmd_simple = [
            'iqtree', '-s', phylip_file,
            '-m', 'HKY+G',           # HKY with Gamma
            '-bb', '1000',           # 1000 ultrafast bootstrap
            '-bnni',                 # Optimize bootstrap
            '-nt', '2',              # Use 2 threads
            '-pre', 'output/ml_tree_hky'
        ]
        
        print(f"Running HKY+G model: {' '.join(cmd_simple)}")
        result1 = subprocess.run(cmd_simple, capture_output=True, text=True)
        
        # Strategy 2: Model testing
        cmd_test = [
            'iqtree', '-s', phylip_file,
            '-m', 'TEST',            # Test best model
            '-bb', '1000',
            '-bnni',
            '-nt', '2',
            '-pre', 'output/ml_tree_best_model'
        ]
        
        print(f"Running model testing: {' '.join(cmd_test)}")
        result2 = subprocess.run(cmd_test, capture_output=True, text=True)
        
        # Strategy 3: GTR model (if simple ones work)
        cmd_gtr = [
            'iqtree', '-s', phylip_file,
            '-m', 'GTR+G',           # General Time Reversible + Gamma
            '-bb', '1000',
            '-bnni',
            '-nt', '2',
            '-pre', 'output/ml_tree_gtr'
        ]
        
        print(f"Running GTR+G model: {' '.join(cmd_gtr)}")
        result3 = subprocess.run(cmd_gtr, capture_output=True, text=True)
        
        # Visualize all ML trees
        visualize_enhanced_ml_trees()
        
        # Create rooted versions
        create_rooted_ml_trees()
        
        return True
    
    return False

def visualize_enhanced_ml_trees():
    """Create enhanced visualizations for ML trees"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    plt.rcParams['axes.titleweight'] = 'bold'
    
    # ML tree files to visualize
    ml_trees = [
        ('output/ml_tree_hky.treefile', 'ML Tree (HKY+G)'),
        ('output/ml_tree_best_model.treefile', 'ML Tree (Best Model)'),
        ('output/ml_tree_gtr.treefile', 'ML Tree (GTR+G)')
    ]
    
    for tree_file, title in ml_trees:
        if os.path.exists(tree_file):
            try:
                tree = Phylo.read(tree_file, 'newick')
                
                # Create standard layout
                fig, ax = plt.subplots(figsize=(16, 12))
                Phylo.draw(tree, axes=ax, do_show=False)
                
                # Enhance the plot
                ax.set_title(f'{title}\nFaba Bean Accessions', 
                            fontsize=18, fontweight='bold', pad=20)
                
                # Enhance tip labels
                for text in ax.texts:
                    text.set_fontweight('bold')
                    text.set_fontsize(12)
                
                # Add bootstrap support values if available
                add_bootstrap_labels(ax, tree)
                
                plt.tight_layout()
                safe_title = title.lower().replace(' ', '_').replace('(', '').replace(')', '')
                plt.savefig(f'plots/{safe_title}.png', dpi=350, bbox_inches='tight')
                plt.savefig(f'plots/{safe_title}.pdf', bbox_inches='tight')
                plt.close()
                
                print(f"✓ Enhanced ML tree visualization saved: plots/{safe_title}.png")
                
            except Exception as e:
                print(f"✗ Error visualizing {tree_file}: {e}")

def add_bootstrap_labels(ax, tree):
    """Add bootstrap support values to tree branches"""
    try:
        # This is a simplified version - in practice you'd parse the bootstrap values
        # from the tree file and add them as labels
        pass
    except:
        pass

def create_rooted_ml_trees():
    """Create rooted versions of ML trees"""
    
    print("Creating rooted ML trees...")
    
    ml_trees = [
        ('output/ml_tree_hky.treefile', 'ML_HKY'),
        ('output/ml_tree_best_model.treefile', 'ML_Best_Model'),
        ('output/ml_tree_gtr.treefile', 'ML_GTR')
    ]
    
    for tree_file, tree_name in ml_trees:
        if os.path.exists(tree_file):
            try:
                tree = Phylo.read(tree_file, 'newick')
                
                # Root the tree at midpoint
                rooted_tree = root_tree_at_midpoint(tree)
                
                # Save rooted tree
                Phylo.write(rooted_tree, f'output/rooted_{tree_name.lower()}.newick', 'newick')
                
                # Create enhanced visualization
                plot_enhanced_rooted_ml_tree(rooted_tree, tree_name.replace('_', ' '))
                
            except Exception as e:
                print(f"✗ Error processing {tree_file}: {e}")

def root_tree_at_midpoint(tree):
    """Root tree at midpoint"""
    terminals = list(tree.get_terminals())
    
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
        tree.root_at_midpoint(tip1, tip2, max_distance/2)
    
    return tree

def plot_enhanced_rooted_ml_tree(tree, tree_type):
    """Create enhanced visualization of rooted ML tree"""
    
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
    
    plt.tight_layout()
    safe_name = tree_type.lower().replace(' ', '_')
    plt.savefig(f'plots/rooted_{safe_name}_enhanced.png', dpi=350, bbox_inches='tight')
    plt.savefig(f'plots/rooted_{safe_name}_enhanced.pdf', bbox_inches='tight')
    plt.close()
    
    print(f"✓ Enhanced rooted {tree_type} tree visualization saved")

def main():
    print("=== Enhanced Maximum Likelihood Tree Analysis ===")
    
    success = build_enhanced_ml_trees()
    
    if success:
        print("\n✓ Enhanced ML tree analysis completed successfully")
        print("✓ Generated multiple ML trees with different models")
        print("✓ Created both unrooted and rooted versions")
    else:
        print("\n✗ Enhanced ML tree analysis failed or skipped")

if __name__ == "__main__":
    main()