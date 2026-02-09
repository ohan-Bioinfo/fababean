# scripts/build_ml_from_phylip.py
import os
import subprocess
import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt

def build_maximum_likelihood_trees():
    """Build Maximum Likelihood trees using IQ-TREE from PHYLIP file"""
    
    print("Building Maximum Likelihood trees...")
    
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
        print("WARNING: IQ-TREE not found. Please install with: conda install -c bioconda iqtree")
        return False
    
    if iqtree_available:
        print("Running IQ-TREE for Maximum Likelihood analysis...")
        
        # Use simpler model for small dataset
        cmd_simple = [
            'iqtree', '-s', phylip_file,
            '-m', 'HKY',           # Hasegawa-Kishino-Yano model (simpler than GTR)
            '-bb', '1000',         # 1000 ultrafast bootstrap replicates
            '-nt', '2',            # Use 2 threads
            '-pre', 'output/ml_tree_simple'
        ]
        
        print(f"Running: {' '.join(cmd_simple)}")
        result = subprocess.run(cmd_simple, capture_output=True, text=True)
        
        if result.returncode == 0:
            print("✓ Simple ML tree analysis completed")
        else:
            print(f"✗ Simple ML tree analysis failed: {result.stderr}")
            # Try without model specification
            cmd_basic = [
                'iqtree', '-s', phylip_file,
                '-bb', '1000',
                '-nt', '2',
                '-pre', 'output/ml_tree_basic'
            ]
            print(f"Trying basic command: {' '.join(cmd_basic)}")
            result = subprocess.run(cmd_basic, capture_output=True, text=True)
            if result.returncode == 0:
                print("✓ Basic ML tree analysis completed")
            else:
                print(f"✗ Basic ML tree analysis failed: {result.stderr}")
        
        # Visualize ML trees
        visualize_ml_trees()
        
        return True
    
    return False

def visualize_ml_trees():
    """Visualize Maximum Likelihood trees"""
    
    plt.rcParams['font.family'] = 'Arial'
    
    # Try to read IQ-TREE output files
    tree_files = [
        ('output/ml_tree_simple.treefile', 'ML Tree (HKY)'),
        ('output/ml_tree_basic.treefile', 'ML Tree (Default Model)')
    ]
    
    for tree_file, title in tree_files:
        if os.path.exists(tree_file):
            try:
                tree = Phylo.read(tree_file, 'newick')
                
                fig, ax = plt.subplots(figsize=(15, 10))
                Phylo.draw(tree, axes=ax, do_show=False)
                
                ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
                
                # Extract filename for saving
                safe_title = title.lower().replace(' ', '_').replace('(', '').replace(')', '')
                plt.tight_layout()
                plt.savefig(f'plots/{safe_title}.png', dpi=300, bbox_inches='tight')
                plt.close()
                
                print(f"✓ ML tree visualization saved: plots/{safe_title}.png")
                
            except Exception as e:
                print(f"✗ Error visualizing {tree_file}: {e}")

def main():
    success = build_maximum_likelihood_trees()
    if success:
        print("\n✓ Maximum Likelihood analysis completed successfully")
    else:
        print("\n✗ Maximum Likelihood analysis failed or skipped")

if __name__ == "__main__":
    main()