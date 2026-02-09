import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from data_loader import SNPDataLoader

def create_summary_dashboard(data):
    """Create a comprehensive summary dashboard"""
    
    fig = plt.figure(figsize=(20, 15))
    fig.suptitle('Faba Bean SNP QC Summary Dashboard', fontsize=20, fontweight='bold', y=0.98)
    
    # Create subplot grid
    gs = fig.add_gridspec(3, 3)
    
    # 1. SNP counts through pipeline
    ax1 = fig.add_subplot(gs[0, 0])
    snp_counts = []
    stages = []
    
    if 'raw_bim' in data:
        snp_counts.append(len(data['raw_bim']))
        stages.append('Raw')
    if 'geno05_bim' in data:
        snp_counts.append(len(data['geno05_bim']))
        stages.append('Geno 0.05')
    if 'maf05_bim' in data:
        snp_counts.append(len(data['maf05_bim']))
        stages.append('MAF 0.05')
    if 'pruned_bim' in data:
        snp_counts.append(len(data['pruned_bim']))
        stages.append('LD Pruned')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    ax1.bar(stages, snp_counts, color=colors[:len(stages)], alpha=0.7)
    ax1.set_title('SNP Counts by Stage', fontweight='bold')
    ax1.tick_params(axis='x', rotation=45)
    ax1.set_ylabel('Number of SNPs')
    
    # 2. Chromosome distribution (raw)
    ax2 = fig.add_subplot(gs[0, 1:])
    if 'raw_bim' in data:
        chr_counts = data['raw_bim']['chr'].value_counts().sort_index()
        ax2.bar(chr_counts.index.astype(str), chr_counts.values, color='skyblue', alpha=0.7)
        ax2.set_title('SNP Distribution - Raw Data', fontweight='bold')
        ax2.set_xlabel('Chromosome')
        ax2.set_ylabel('SNP Count')
    
    # 3. Sample missingness
    ax3 = fig.add_subplot(gs[1, 0])
    if 'sample_missingness' in data:
        ax3.hist(data['sample_missingness']['F_MISS'], bins=30, 
                alpha=0.7, color='coral')
        ax3.set_title('Sample Missingness', fontweight='bold')
        ax3.set_xlabel('Missing Rate')
        ax3.set_ylabel('Count')
    
    # 4. SNP missingness
    ax4 = fig.add_subplot(gs[1, 1])
    if 'snp_missingness' in data:
        ax4.hist(data['snp_missingness']['F_MISS'], bins=30, 
                alpha=0.7, color='lightgreen')
        ax4.set_title('SNP Missingness', fontweight='bold')
        ax4.set_xlabel('Missing Rate')
        ax4.set_ylabel('Count')
    
    # 5. Heterozygosity
    ax5 = fig.add_subplot(gs[1, 2])
    if 'heterozygosity' in data and 'O(HOM)' in data['heterozygosity'].columns and 'N(NM)' in data['heterozygosity'].columns:
        het_data = data['heterozygosity']
        het_rates = (het_data['N(NM)'] - het_data['O(HOM)']) / het_data['N(NM)']
        ax5.hist(het_rates, bins=30, alpha=0.7, color='mediumpurple')
        ax5.set_title('Heterozygosity Rates', fontweight='bold')
        ax5.set_xlabel('Heterozygosity Rate')
        ax5.set_ylabel('Count')
    
    # 6. Filtering efficiency
    ax6 = fig.add_subplot(gs[2, :])
    if len(snp_counts) > 0 and 'Raw' in stages:
        raw_index = stages.index('Raw')
        retention = [100]
        for i, stage in enumerate(stages):
            if i != raw_index:
                retention.append((snp_counts[i]/snp_counts[raw_index])*100)
        
        ax6.plot([stages[raw_index]] + [s for i, s in enumerate(stages) if i != raw_index], 
                retention, marker='o', linewidth=3, markersize=8)
        ax6.set_title('SNP Retention Percentage', fontweight='bold')
        ax6.set_ylabel('Retention (%)')
        ax6.set_xlabel('Filtering Stage')
        ax6.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_QC_Dashboard.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    create_summary_dashboard(data)

