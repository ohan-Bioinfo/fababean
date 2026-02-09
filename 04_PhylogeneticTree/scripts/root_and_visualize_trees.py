# PhylogeneticTree/scripts/root_and_visualize_trees.py
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import matplotlib.pyplot as plt
import numpy as np

def root_trees_with_midpoint():
    """Apply midpoint rooting to all trees"""
    
    print("Applying midpoint rooting to trees...")
    
    # List of tree files to root
    tree_files = [
        ('output/nj_tree_euclidean.newick', 'NJ Euclidean'),
        ('output/nj_tree_hamming.newick', 'NJ Hamming'),
        ('output/ml_tree_basic.treefile', 'ML Basic'),
        ('output/ml_tree_comprehensive.treefile', 'ML Comprehensive')
    ]
    
    rooted_trees = []
    
    for tree_file, tree_name in tree_files:
        try:
            # Read tree
            tree = Phylo.read(tree_file, 'newick')
            
            # Apply midpoint rooting
            rooted_tree = root_tree_at_midpoint(tree)
            
            # Save rooted tree
            rooted_filename = f'output/rooted_{tree_name.lower().replace(" ", "_")}.newick'
            Phylo.write(rooted_tree, rooted_filename, 'newick')
            
            rooted_trees.append((rooted_tree, tree_name, rooted_filename))
            
            print(f"✓ Rooted {tree_name} tree")
            
        except Exception as e:
            print(f"✗ Could not root {tree_name}: {e}")
    
    # Create comprehensive visualization
    create_comprehensive_visualization(rooted_trees)
    
    return rooted_trees

def root_tree_at_midpoint(tree):
    """Root tree at midpoint (longest path midpoint)"""
    
    # Find the two most distant tips
    max_distance = 0
    tip1, tip2 = None, None
    
    for tip1_candidate in tree.get_terminals():
        for tip2_candidate in tree.get_terminals():
            if tip1_candidate != tip2_candidate:
                distance = tree.distance(tip1_candidate, tip2_candidate)
                if distance > max_distance:
                    max_distance = distance
                    tip1, tip2 = tip1_candidate, tip2_candidate
    
    if tip1 and tip2:
        # Find midpoint
        midpoint = max_distance / 2
        
        # Root the tree at the midpoint
        tree.root_at_midpoint(tip1, tip2, midpoint)
    
    return tree

def create_comprehensive_visualization(rooted_trees):
    """Create comprehensive visualization of all rooted trees"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 10
    
    n_trees = len(rooted_trees)
    if n_trees == 0:
        print("No rooted trees to visualize")
        return
    
    # Create subplot grid
    if n_trees <= 2:
        fig, axes = plt.subplots(1, n_trees, figsize=(6*n_trees, 8))
        if n_trees == 1:
            axes = [axes]
    else:
        rows = (n_trees + 1) // 2
        fig, axes = plt.subplots(rows, 2, figsize=(16, 4*rows))
        axes = axes.flatten()
    
    # Plot each tree
    for i, (tree, tree_name, filename) in enumerate(rooted_trees):
        if i < len(axes):
            ax = axes[i]
            Phylo.draw(tree, axes=ax, do_show=False)
            ax.set_title(f'Rooted {tree_name} Tree\n(Midpoint Rooting)', 
                        fontsize=12, fontweight='bold', pad=20)
    
    # Hide unused subplots
    for i in range(n_trees, len(axes)):
        axes[i].set_visible(False)
    
    plt.tight_layout()
    plt.savefig('plots/all_rooted_trees_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('plots/all_rooted_trees_comparison.pdf', bbox_inches='tight')
    
    # Create individual high-quality plots for each tree
    for tree, tree_name, filename in rooted_trees:
        fig, ax = plt.subplots(figsize=(12, 8))
        Phylo.draw(tree, axes=ax, do_show=False)
        
        # Customize appearance
        ax.set_title(f'Rooted {tree_name} Phylogenetic Tree\nFaba Bean Accessions (Midpoint Rooting)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add summary text
        n_tips = len(tree.get_terminals())
        ax.text(0.02, 0.98, f'Total tips: {n_tips}\nRooting: Midpoint', 
               transform=ax.transAxes, fontsize=10, 
               verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        safe_name = tree_name.lower().replace(' ', '_')
        plt.savefig(f'plots/rooted_{safe_name}_detailed.png', dpi=300, bbox_inches='tight')
        plt.close()
    
    print("Comprehensive tree visualizations saved")

def generate_phylogenetic_summary(rooted_trees):
    """Generate summary statistics for phylogenetic analysis"""
    
    print("\n=== PHYLOGENETIC ANALYSIS SUMMARY ===")
    
    for tree, tree_name, filename in rooted_trees:
        print(f"\n--- {tree_name} Tree ---")
        print(f"Number of tips: {len(tree.get_terminals())}")
        print(f"Number of internal nodes: {len(tree.get_nonterminals())}")
        print(f"Tree depth: {tree.depths().max():.4f}")
        print(f"Tree file: {filename}")
    
    # Save summary to file
    with open('output/phylogenetic_analysis_summary.txt', 'w') as f:
        f.write("Faba Bean Phylogenetic Analysis Summary\n")
        f.write("======================================\n\n")
        
        f.write("Dataset: 150 fingerprint SNPs across 21 accessions\n\n")
        
        f.write("Methods Applied:\n")
        f.write("- Neighbour-Joining (Euclidean distance)\n")
        f.write("- Neighbour-Joining (Hamming distance)\n") 
        f.write("- Maximum Likelihood (GTR+G model)\n")
        f.write("- Maximum Likelihood (Model testing)\n\n")
        
        f.write("Rooting: Midpoint rooting (no outgroup available)\n\n")
        
        f.write("Output Files:\n")
        for _, tree_name, filename in rooted_trees:
            f.write(f"- {tree_name}: {filename}\n")

if __name__ == "__main__":
    rooted_trees = root_trees_with_midpoint()
    generate_phylogenetic_summary(rooted_trees)
    
    print("\n=== PHYLOGENETIC ANALYSIS COMPLETE ===")
    print("✓ Unrooted trees created (NJ + ML methods)")
    print("✓ Rooted trees created (midpoint rooting)") 
    print("✓ Comprehensive visualizations generated")
    print("✓ Summary statistics calculated")
    print("✓ All output files saved in PhylogeneticTree/ directory")