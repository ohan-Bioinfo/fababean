import matplotlib.pyplot as plt
import pandas as pd
from data_loader import SNPDataLoader

def plot_filtering_pipeline(data):
    """Plot SNP counts through the filtering pipeline"""
    
    snp_counts = {}
    if 'raw_bim' in data:
        snp_counts['Raw'] = len(data['raw_bim'])
    if 'geno05_bim' in data:
        snp_counts['Geno 0.05'] = len(data['geno05_bim'])
    if 'maf05_bim' in data:
        snp_counts['MAF 0.05'] = len(data['maf05_bim'])
    if 'pruned_bim' in data:
        snp_counts['LD Pruned'] = len(data['pruned_bim'])
    
    if not snp_counts:
        print("No data available for filtering pipeline")
        return
    
    # Plot 1: Bar plot of SNP counts
    fig1, ax1 = plt.subplots(figsize=(12, 6))
    stages = list(snp_counts.keys())
    counts = list(snp_counts.values())
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    bars = ax1.bar(stages, counts, color=colors[:len(stages)], alpha=0.7)
    
    # Add value labels on bars
    for bar, count in zip(bars, counts):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + max(counts)*0.01,
               f'{count:,}', ha='center', va='bottom', fontweight='bold')
    
    ax1.set_title('SNP Counts Through Filtering Pipeline', fontsize=16, fontweight='bold')
    ax1.set_ylabel('Number of SNPs')
    ax1.set_xlabel('Filtering Stage')
    ax1.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_Filtering_Pipeline.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Plot 2: Retention percentage
    if 'Raw' in snp_counts:
        fig2, ax2 = plt.subplots(figsize=(12, 6))
        retention_pct = [100]
        for stage in stages[1:]:
            retention_pct.append((snp_counts[stage]/snp_counts['Raw'])*100)
        
        ax2.plot(stages, retention_pct, marker='o', linewidth=2, markersize=8, 
                color='purple', markerfacecolor='red')
        ax2.set_title('SNP Retention Percentage Through Filtering Pipeline', 
                     fontsize=16, fontweight='bold')
        ax2.set_ylabel('Retention (%)')
        ax2.set_xlabel('Filtering Stage')
        ax2.grid(True, alpha=0.3)
        
        # Add percentage labels
        for i, pct in enumerate(retention_pct):
            ax2.annotate(f'{pct:.1f}%', (stages[i], pct), 
                        textcoords="offset points", xytext=(0,10), 
                        ha='center', fontweight='bold')
        
        plt.tight_layout()
        plt.savefig('Faba_Bean_Retention_Percentage.png', dpi=300, bbox_inches='tight')
        plt.show()

if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    plot_filtering_pipeline(data)
