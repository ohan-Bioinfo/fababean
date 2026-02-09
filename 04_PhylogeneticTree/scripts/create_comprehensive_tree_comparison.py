# scripts/create_comprehensive_tree_comparison.py
from Bio import Phylo
import matplotlib.pyplot as plt
import os
import numpy as np

def create_comprehensive_comparison():
    """Create comprehensive comparison of all phylogenetic trees"""
    
    print("Creating comprehensive tree comparison...")
    
    # List all tree files
    tree_files = []
    
    # NJ trees
    if os.path.exists('output/nj_tree_unrooted.newick'):
        tree_files.append(('output/nj_tree_unrooted.newick', 'NJ Unrooted'))
    if os.path.exists('output/nj_tree_rooted.newick'):
        tree_files.append(('output/nj_tree_rooted.newick', 'NJ Rooted'))
    
    # ML trees
    ml_trees = [
        ('output/ml_tree_hky.treefile', 'ML HKY'),
        ('output/ml_tree_best_model.treefile', 'ML Best Model'),
        ('output/ml_tree_gtr.treefile', 'ML GTR'),
        ('output/rooted_ml_hky.newick', 'ML HKY Rooted'),
        ('output/rooted_ml_best_model.newick', 'ML Best Model Rooted'),
        ('output/rooted_ml_gtr.newick', 'ML GTR Rooted')
    ]
    
    for tree_file, tree_name in ml_trees:
        if os.path.exists(tree_file):
            tree_files.append((tree_file, tree_name))
    
    print(f"Found {len(tree_files)} tree files for comparison")
    
    # Create comparison figure
    if tree_files:
        create_tree_comparison_figure(tree_files)
        create_individual_tree_plots(tree_files)
        generate_tree_comparison_summary(tree_files)
    
    return tree_files

def create_tree_comparison_figure(tree_files):
    """Create a figure comparing all trees"""
    
    n_trees = len(tree_files)
    if n_trees == 0:
        return
    
    # Calculate grid size
    cols = min(2, n_trees)
    rows = (n_trees + cols - 1) // cols
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    
    fig, axes = plt.subplots(rows, cols, figsize=(6*cols, 5*rows))
    
    if n_trees == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for i, (tree_file, tree_name) in enumerate(tree_files):
        if i < len(axes):
            try:
                tree = Phylo.read(tree_file, 'newick')
                ax = axes[i]
                
                Phylo.draw(tree, axes=ax, do_show=False)
                ax.set_title(tree_name, fontsize=14, fontweight='bold', pad=10)
                
                # Enhance tip labels
                for text in ax.texts:
                    text.set_fontweight('bold')
                    text.set_fontsize(8)
                
            except Exception as e:
                print(f"✗ Error plotting {tree_name}: {e}")
                axes[i].set_title(f"{tree_name}\n(Error)", color='red')
                axes[i].text(0.5, 0.5, "Plotting Error", 
                           ha='center', va='center', transform=axes[i].transAxes)
    
    # Hide unused subplots
    for i in range(n_trees, len(axes)):
        axes[i].set_visible(False)
    
    plt.suptitle('Comprehensive Phylogenetic Tree Comparison\nFaba Bean Accessions', 
                fontsize=16, fontweight='bold', y=0.95)
    plt.tight_layout()
    plt.savefig('plots/comprehensive_tree_comparison.png', dpi=300, bbox_inches='tight')
    plt.savefig('plots/comprehensive_tree_comparison.pdf', bbox_inches='tight')
    plt.close()
    
    print("✓ Comprehensive tree comparison figure saved")

def create_individual_tree_plots(tree_files):
    """Create individual high-quality plots for each tree"""
    
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.weight'] = 'bold'
    
    for tree_file, tree_name in tree_files:
        try:
            tree = Phylo.read(tree_file, 'newick')
            
            fig, ax = plt.subplots(figsize=(12, 8))
            Phylo.draw(tree, axes=ax, do_show=False)
            
            # Enhanced title and labels
            title_suffix = " (Rooted)" if "Rooted" in tree_name else " (Unrooted)"
            ax.set_title(f'{tree_name}{title_suffix}\nFaba Bean Accessions', 
                        fontsize=16, fontweight='bold', pad=20)
            
            # Enhanced tip labels
            for text in ax.texts:
                text.set_fontweight('bold')
                text.set_fontsize(10)
            
            # Add tree statistics
            n_tips = len(tree.get_terminals())
            tree_length = sum(tree.depths().values())
            
            ax.text(0.02, 0.98, f'Samples: {n_tips}\nTree length: {tree_length:.2f}', 
                   transform=ax.transAxes, fontsize=9, fontweight='bold',
                   verticalalignment='top', 
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
            
            plt.tight_layout()
            safe_name = tree_name.lower().replace(' ', '_')
            plt.savefig(f'plots/{safe_name}_detailed.png', dpi=350, bbox_inches='tight')
            plt.savefig(f'plots/{safe_name}_detailed.pdf', bbox_inches='tight')
            plt.close()
            
            print(f"✓ Individual tree plot saved: plots/{safe_name}_detailed.png")
            
        except Exception as e:
            print(f"✗ Error creating individual plot for {tree_name}: {e}")

def generate_tree_comparison_summary(tree_files):
    """Generate summary of all trees"""
    
    print("\n" + "="*60)
    print("COMPREHENSIVE PHYLOGENETIC ANALYSIS SUMMARY")
    print("="*60)
    
    summary_data = []
    
    for tree_file, tree_name in tree_files:
        try:
            tree = Phylo.read(tree_file, 'newick')
            n_tips = len(tree.get_terminals())
            tree_length = sum(tree.depths().values())
            
            summary_data.append({
                'Tree': tree_name,
                'Samples': n_tips,
                'Tree_Length': f"{tree_length:.4f}",
                'Method': 'NJ' if 'NJ' in tree_name else 'ML',
                'Rooted': 'Yes' if 'Rooted' in tree_name else 'No'
            })
            
            print(f"\n--- {tree_name} ---")
            print(f"  Samples: {n_tips}")
            print(f"  Tree length: {tree_length:.4f}")
            print(f"  Method: {'NJ' if 'NJ' in tree_name else 'ML'}")
            print(f"  Rooted: {'Yes' if 'Rooted' in tree_name else 'No'}")
            
        except Exception as e:
            print(f"\n--- {tree_name} ---")
            print(f"  ERROR: {e}")
    
    # Save summary to CSV
    if summary_data:
        import pandas as pd
        df_summary = pd.DataFrame(summary_data)
        df_summary.to_csv('output/phylogenetic_trees_summary.csv', index=False)
        print(f"\n✓ Tree summary saved: output/phylogenetic_trees_summary.csv")
    
    print(f"\nTotal trees generated: {len(tree_files)}")
    print("✓ Analysis complete!")

def main():
    print("=== Comprehensive Phylogenetic Tree Comparison ===")
    
    tree_files = create_comprehensive_comparison()
    
    if tree_files:
        print(f"\n✓ Successfully processed {len(tree_files)} trees")
        print("✓ Comprehensive comparison completed")
    else:
        print("\n✗ No trees found for comparison")

if __name__ == "__main__":
    main()