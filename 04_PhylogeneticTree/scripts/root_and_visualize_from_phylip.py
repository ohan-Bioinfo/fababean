# scripts/root_and_visualize_from_phylip.py
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
import matplotlib.pyplot as plt
import numpy as np

def root_trees_with_midpoint():
    """Apply midpoint rooting to all trees"""
    
    print("Applying midpoint rooting to trees...")
    
    # List of tree files to root
    tree_files = [
        ('output/nj_tree_unrooted.newick', 'NJ'),
        ('output/ml_tree_basic.treefile', 'ML Basic'),
        ('output/ml_tree_model_test.treefile', 'ML Model Test')
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
    if rooted_trees:
        create_comprehensive_visualization(rooted_trees)
        generate_analysis_summary(rooted_trees)
    
    return rooted_trees

def root_tree_at_midpoint(tree):
    """Root tree at midpoint (longest path midpoint)"""
    
    # Find the two most distant tips
    max_distance = 0
    tip1, tip2 = None, None
    
    terminals = list(tree.get_terminals())
    
    for tip1_candidate in terminals:
        for tip2_candidate in terminals:
            if tip1_candidate != tip2_candidate:
                distance = tree.distance(tip1_candidate, tip2_candidate)
                if distance > max_distance:
                    max_distance = distance
                    tip1, tip2 = tip1_candidate, tip2_candidate
    
    if tip1 and tip2:
        # Root at midpoint
        tree.root_at_midpoint(tip1, tip2, max_distance/2)
    
    return tree

def create_comprehensive_visualization(rooted_trees):
    """Create comprehensive visualization of all rooted trees"""
    
    plt.rcParams['font.family'] = 'Arial'
    
    n_trees = len(rooted_trees)
    
    # Create individual plots for each tree
    for tree, tree_name, filename in rooted_trees:
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Draw tree with better formatting
        Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: x.name)
        
        # Customize appearance
        ax.set_title(f'Rooted {tree_name} Phylogenetic Tree\nFaba Bean Accessions (Midpoint Rooting)', 
                    fontsize=14, fontweight='bold', pad=20)
        
        # Add tree statistics
        n_tips = len(tree.get_terminals())
        tree_length = sum(tree.depths().values())
        
        ax.text(0.02, 0.98, f'Samples: {n_tips}\nTree length: {tree_length:.2f}', 
               transform=ax.transAxes, fontsize=10, 
               verticalalignment='top', 
               bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
        
        plt.tight_layout()
        safe_name = tree_name.lower().replace(' ', '_')
        plt.savefig(f'plots/rooted_{safe_name}_detailed.png', dpi=300, bbox_inches='tight')
        plt.savefig(f'plots/rooted_{safe_name}_detailed.pdf', bbox_inches='tight')
        plt.close()
    
    print("✓ Rooted tree visualizations saved")

def generate_analysis_summary(rooted_trees):
    """Generate summary of phylogenetic analysis"""
    
    print("\n" + "="*50)
    print("PHYLOGENETIC ANALYSIS SUMMARY")
    print("="*50)
    
    for tree, tree_name, filename in rooted_trees:
        print(f"\n--- {tree_name} Tree ---")
        print(f"Number of tips: {len(tree.get_terminals())}")
        print(f"Tree length: {sum(tree.depths().values()):.4f}")
        print(f"Tree file: {filename}")
    
    # Save summary to file
    with open('output/phylogenetic_analysis_summary.txt', 'w') as f:
        f.write("Faba Bean Phylogenetic Analysis Summary\n")
        f.write("======================================\n\n")
        
        f.write("Dataset: 150 fingerprint SNPs across 21 accessions\n")
        f.write("Input format: VCF -> PHYLIP (via vcf2phylip.py)\n\n")
        
        f.write("Methods Applied:\n")
        f.write("- Neighbour-Joining (Hamming distance)\n")
        f.write("- Maximum Likelihood (GTR+G model)\n")
        f.write("- Maximum Likelihood (Model testing)\n\n")
        
        f.write("Rooting: Midpoint rooting\n\n")
        
        f.write("Output Files:\n")
        for _, tree_name, filename in rooted_trees:
            f.write(f"- {tree_name}: {filename}\n")

def main():
    print("=== Final Phylogenetic Tree Processing ===")
    
    rooted_trees = root_trees_with_midpoint()
    
    if rooted_trees:
        print(f"\n✓ Successfully processed {len(rooted_trees)} trees")
        print("✓ Rooted trees saved in output/ directory")
        print("✓ Visualizations saved in plots/ directory")
        print("✓ Analysis summary generated")
    else:
        print("\n✗ No trees were successfully processed")
    
    print("\n=== Phylogenetic Analysis Complete ===")

if __name__ == "__main__":
    main()