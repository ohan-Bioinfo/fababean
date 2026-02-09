import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from data_loader import SNPDataLoader

def plot_heterozygosity(data):
    """Plot heterozygosity analysis with proper column handling"""
    
    if 'heterozygosity' not in data:
        print("Heterozygosity data not available")
        return
    
    het_data = data['heterozygosity']
    
    # Print available columns for debugging
    print(f"Available columns in heterozygosity data: {het_data.columns.tolist()}")
    
    # Based on your output, columns are: ['FID', 'IID', 'O(HOM)', 'E(HOM)', 'N(NM)', 'F']
    # N(NM) is the number of non-missing genotypes, which we can use as total sites
    if 'O(HOM)' in het_data.columns and 'N(NM)' in het_data.columns:
        het_data['HET_RATE'] = (het_data['N(NM)'] - het_data['O(HOM)']) / het_data['N(NM)']
        print("✓ Heterozygosity rate calculated using O(HOM) and N(NM)")
    else:
        print("Required columns for heterozygosity calculation not found")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    fig.suptitle('Sample Heterozygosity Analysis', fontsize=16, fontweight='bold')
    
    # 1. Heterozygosity rate distribution
    axes[0,0].hist(het_data['HET_RATE'], bins=20, alpha=0.7, 
                  color='mediumpurple', edgecolor='black')
    axes[0,0].set_title('A. Heterozygosity Rate Distribution', fontweight='bold')
    axes[0,0].set_xlabel('Heterozygosity Rate')
    axes[0,0].set_ylabel('Number of Samples')
    axes[0,0].grid(alpha=0.3)
    
    # 2. Heterozygosity vs Homozygosity
    scatter = axes[0,1].scatter(het_data['O(HOM)'], het_data['HET_RATE'], 
                               alpha=0.7, c=het_data['N(NM)'], 
                               cmap='viridis', s=60)
    axes[0,1].set_title('B. Heterozygosity vs Homozygous Sites', fontweight='bold')
    axes[0,1].set_xlabel('Number of Homozygous Sites [O(HOM)]')
    axes[0,1].set_ylabel('Heterozygosity Rate')
    axes[0,1].grid(alpha=0.3)
    plt.colorbar(scatter, ax=axes[0,1], label='Non-missing Sites [N(NM)]')
    
    # 3. Inbreeding coefficient distribution
    if 'F' in het_data.columns:
        axes[1,0].hist(het_data['F'], bins=20, alpha=0.7, 
                      color='orange', edgecolor='black')
        axes[1,0].set_title('C. Inbreeding Coefficient (F) Distribution', fontweight='bold')
        axes[1,0].set_xlabel('Inbreeding Coefficient (F)')
        axes[1,0].set_ylabel('Number of Samples')
        axes[1,0].grid(alpha=0.3)
    
    # 4. Expected vs Observed homozygosity
    if 'E(HOM)' in het_data.columns:
        axes[1,1].scatter(het_data['E(HOM)'], het_data['O(HOM)'], 
                         alpha=0.7, color='red', s=60)
        axes[1,1].plot([het_data['E(HOM)'].min(), het_data['E(HOM)'].max()], 
                      [het_data['E(HOM)'].min(), het_data['E(HOM)'].max()], 
                      'k--', alpha=0.5, label='y=x')
        axes[1,1].set_title('D. Observed vs Expected Homozygosity', fontweight='bold')
        axes[1,1].set_xlabel('Expected Homozygous Sites [E(HOM)]')
        axes[1,1].set_ylabel('Observed Homozygous Sites [O(HOM)]')
        axes[1,1].legend()
        axes[1,1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_Heterozygosity_Analysis.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print heterozygosity statistics
    print(f"\nHeterozygosity Statistics:")
    print(f"Mean heterozygosity rate: {het_data['HET_RATE'].mean():.4f}")
    print(f"Std heterozygosity rate: {het_data['HET_RATE'].std():.4f}")
    print(f"Min heterozygosity rate: {het_data['HET_RATE'].min():.4f}")
    print(f"Max heterozygosity rate: {het_data['HET_RATE'].max():.4f}")

if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    plot_heterozygosity(data)
