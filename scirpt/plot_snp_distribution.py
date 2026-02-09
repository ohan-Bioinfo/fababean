import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from data_loader import SNPDataLoader

def plot_snp_distribution_main_chromosomes(data):
    """Plot SNP distribution across the 7 main chromosomes for different filtering stages"""
    
    if 'raw_bim' not in data:
        print("Error: BIM files not loaded")
        return
    
    # Define the main chromosomes based on your data
    main_chromosomes = ['chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']
    
    # Filter for main chromosomes
    def filter_main_chromosomes(bim_df):
        bim_df = bim_df.copy()
        bim_df['chr_str'] = bim_df['chr'].astype(str)
        return bim_df[bim_df['chr_str'].isin(main_chromosomes)]
    
    # Apply filtering to all BIM files
    raw_filtered = filter_main_chromosomes(data['raw_bim'])
    geno05_filtered = filter_main_chromosomes(data['geno05_bim']) if 'geno05_bim' in data else None
    maf05_filtered = filter_main_chromosomes(data['maf05_bim']) if 'maf05_bim' in data else None
    pruned_filtered = filter_main_chromosomes(data['pruned_bim']) if 'pruned_bim' in data else None
    
    # Create the plot
    fig, axes = plt.subplots(2, 2, figsize=(20, 12))
    fig.suptitle('SNP Distribution Across Main Chromosomes at Different Filtering Stages', 
                fontsize=16, fontweight='bold')
    
    # Plot 1: Raw data
    chr_counts_raw = raw_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    bars_raw = axes[0,0].bar(range(len(chr_counts_raw)), chr_counts_raw.values, 
                           color='skyblue', alpha=0.7, edgecolor='black')
    axes[0,0].set_title('A. Raw SNPs', fontweight='bold')
    axes[0,0].set_xlabel('Chromosome')
    axes[0,0].set_ylabel('Number of SNPs')
    axes[0,0].set_xticks(range(len(chr_counts_raw)))
    axes[0,0].set_xticklabels(chr_counts_raw.index, rotation=45)
    
    # Add numbers on top of bars for raw data
    for bar, count in zip(bars_raw, chr_counts_raw.values):
        height = bar.get_height()
        axes[0,0].text(bar.get_x() + bar.get_width()/2., height + max(chr_counts_raw.values)*0.01,
                      f'{int(count):,}', ha='center', va='bottom', fontweight='bold', fontsize=10)
    
    # Plot 2: After genotype filtering
    if geno05_filtered is not None:
        chr_counts_geno05 = geno05_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
        bars_geno05 = axes[0,1].bar(range(len(chr_counts_geno05)), chr_counts_geno05.values, 
                                  color='lightcoral', alpha=0.7, edgecolor='black')
        axes[0,1].set_title('B. After Genotype Call Rate Filter (--geno 0.05)', fontweight='bold')
        axes[0,1].set_xlabel('Chromosome')
        axes[0,1].set_ylabel('Number of SNPs')
        axes[0,1].set_xticks(range(len(chr_counts_geno05)))
        axes[0,1].set_xticklabels(chr_counts_geno05.index, rotation=45)
        
        # Add numbers on top of bars
        for bar, count in zip(bars_geno05, chr_counts_geno05.values):
            height = bar.get_height()
            axes[0,1].text(bar.get_x() + bar.get_width()/2., height + max(chr_counts_geno05.values)*0.01,
                          f'{int(count):,}', ha='center', va='bottom', fontweight='bold', fontsize=10)
    
    # Plot 3: After MAF filtering
    if maf05_filtered is not None:
        chr_counts_maf05 = maf05_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
        bars_maf05 = axes[1,0].bar(range(len(chr_counts_maf05)), chr_counts_maf05.values, 
                                 color='lightgreen', alpha=0.7, edgecolor='black')
        axes[1,0].set_title('C. After MAF Filter (MAF > 0.05)', fontweight='bold')
        axes[1,0].set_xlabel('Chromosome')
        axes[1,0].set_ylabel('Number of SNPs')
        axes[1,0].set_xticks(range(len(chr_counts_maf05)))
        axes[1,0].set_xticklabels(chr_counts_maf05.index, rotation=45)
        
        # Add numbers on top of bars
        for bar, count in zip(bars_maf05, chr_counts_maf05.values):
            height = bar.get_height()
            axes[1,0].text(bar.get_x() + bar.get_width()/2., height + max(chr_counts_maf05.values)*0.01,
                          f'{int(count):,}', ha='center', va='bottom', fontweight='bold', fontsize=10)
    
    # Plot 4: After LD pruning
    if pruned_filtered is not None:
        chr_counts_pruned = pruned_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
        bars_pruned = axes[1,1].bar(range(len(chr_counts_pruned)), chr_counts_pruned.values, 
                                  color='gold', alpha=0.7, edgecolor='black')
        axes[1,1].set_title('D. After LD Pruning', fontweight='bold')
        axes[1,1].set_xlabel('Chromosome')
        axes[1,1].set_ylabel('Number of SNPs')
        axes[1,1].set_xticks(range(len(chr_counts_pruned)))
        axes[1,1].set_xticklabels(chr_counts_pruned.index, rotation=45)
        
        # Add numbers on top of bars
        for bar, count in zip(bars_pruned, chr_counts_pruned.values):
            height = bar.get_height()
            axes[1,1].text(bar.get_x() + bar.get_width()/2., height + max(chr_counts_pruned.values)*0.01,
                          f'{int(count):,}', ha='center', va='bottom', fontweight='bold', fontsize=10)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_SNP_Distribution_Main_Chromosomes.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print summary statistics
    print("\n=== SNP Counts Summary (Main Chromosomes) ===")
    print(f"Raw SNPs: {len(raw_filtered):,}")
    if geno05_filtered is not None:
        print(f"After genotype filter: {len(geno05_filtered):,}")
    if maf05_filtered is not None:
        print(f"After MAF filter: {len(maf05_filtered):,}")
    if pruned_filtered is not None:
        print(f"After LD pruning: {len(pruned_filtered):,}")
    
    # Print detailed breakdown
    print("\n=== Detailed Breakdown by Chromosome ===")
    for chrom in main_chromosomes:
        raw_count = len(raw_filtered[raw_filtered['chr_str'] == chrom])
        pruned_count = len(pruned_filtered[pruned_filtered['chr_str'] == chrom]) if pruned_filtered is not None else 0
        retention = (pruned_count / raw_count * 100) if raw_count > 0 else 0
        print(f"{chrom}: Raw={raw_count:,}, Pruned={pruned_count:,}, Retention={retention:.1f}%")

def plot_chromosome_comparison(data):
    """Compare SNP counts across main chromosomes for different stages"""
    if 'raw_bim' not in data:
        return
    
    # Define main chromosomes
    main_chromosomes = ['chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']
    
    # Filter for main chromosomes
    def filter_main_chromosomes(bim_df):
        bim_df = bim_df.copy()
        bim_df['chr_str'] = bim_df['chr'].astype(str)
        return bim_df[bim_df['chr_str'].isin(main_chromosomes)]
    
    # Get chromosome counts for each stage
    chr_counts_data = {}
    raw_filtered = filter_main_chromosomes(data['raw_bim'])
    chr_counts_data['Raw'] = raw_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    
    if 'geno05_bim' in data:
        geno05_filtered = filter_main_chromosomes(data['geno05_bim'])
        chr_counts_data['Geno_0.05'] = geno05_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    
    if 'maf05_bim' in data:
        maf05_filtered = filter_main_chromosomes(data['maf05_bim'])
        chr_counts_data['MAF_0.05'] = maf05_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    
    if 'pruned_bim' in data:
        pruned_filtered = filter_main_chromosomes(data['pruned_bim'])
        chr_counts_data['LD_Pruned'] = pruned_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    
    # Create DataFrame with consistent chromosome order
    chr_counts = pd.DataFrame(chr_counts_data)
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x = np.arange(len(chr_counts.index))
    width = 0.2
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    for i, (stage, color) in enumerate(zip(chr_counts.columns, colors)):
        if i < len(chr_counts.columns):
            offset = width * (i - (len(chr_counts.columns)-1)/2)
            bars = ax.bar(x + offset, chr_counts[stage], width, label=stage, alpha=0.7, color=color, edgecolor='black')
            
            # Add numbers on top of bars
            for bar, count in zip(bars, chr_counts[stage]):
                if count > 0:  # Only add text if count > 0
                    height = bar.get_height()
                    ax.text(bar.get_x() + bar.get_width()/2., height + max(chr_counts[stage])*0.01,
                           f'{int(count):,}', ha='center', va='bottom', fontweight='bold', fontsize=8)
    
    ax.set_title('SNP Counts by Chromosome Across Filtering Stages', 
                fontsize=16, fontweight='bold')
    ax.set_xlabel('Chromosome')
    ax.set_ylabel('Number of SNPs')
    ax.set_xticks(x)
    ax.set_xticklabels(chr_counts.index)
    ax.legend()
    ax.grid(axis='y', alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_Chromosome_Comparison_Main.png', dpi=300, bbox_inches='tight')
    plt.show()

def plot_snp_retention_by_chromosome(data):
    """Plot SNP retention percentage by chromosome"""
    if 'raw_bim' not in data or 'pruned_bim' not in data:
        print("Raw or pruned data not available for retention plot")
        return
    
    # Define main chromosomes
    main_chromosomes = ['chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']
    
    # Filter for main chromosomes
    def filter_main_chromosomes(bim_df):
        bim_df = bim_df.copy()
        bim_df['chr_str'] = bim_df['chr'].astype(str)
        return bim_df[bim_df['chr_str'].isin(main_chromosomes)]
    
    raw_filtered = filter_main_chromosomes(data['raw_bim'])
    pruned_filtered = filter_main_chromosomes(data['pruned_bim'])
    
    # Get counts by chromosome
    raw_counts = raw_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    pruned_counts = pruned_filtered['chr_str'].value_counts().reindex(main_chromosomes).fillna(0)
    
    # Calculate retention percentage
    retention_pct = (pruned_counts / raw_counts * 100).fillna(0)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 6))
    
    # Plot 1: Absolute counts
    x = np.arange(len(raw_counts.index))
    width = 0.35
    
    bars1 = ax1.bar(x - width/2, raw_counts.values, width, label='Raw SNPs', color='skyblue', alpha=0.7, edgecolor='black')
    bars2 = ax1.bar(x + width/2, pruned_counts.values, width, label='LD Pruned SNPs', color='gold', alpha=0.7, edgecolor='black')
    
    ax1.set_title('SNP Counts: Raw vs LD Pruned by Chromosome', fontweight='bold')
    ax1.set_xlabel('Chromosome')
    ax1.set_ylabel('Number of SNPs')
    ax1.set_xticks(x)
    ax1.set_xticklabels(raw_counts.index)
    ax1.legend()
    ax1.grid(axis='y', alpha=0.3)
    
    # Add numbers on bars
    for bars in [bars1, bars2]:
        for bar in bars:
            height = bar.get_height()
            if height > 0:
                ax1.text(bar.get_x() + bar.get_width()/2., height + max(raw_counts.values)*0.01,
                        f'{int(height):,}', ha='center', va='bottom', fontweight='bold', fontsize=8)
    
    # Plot 2: Retention percentage
    bars_retention = ax2.bar(range(len(retention_pct)), retention_pct.values, 
                           color='lightgreen', alpha=0.7, edgecolor='black')
    ax2.set_title('SNP Retention Percentage After LD Pruning', fontweight='bold')
    ax2.set_xlabel('Chromosome')
    ax2.set_ylabel('Retention (%)')
    ax2.set_xticks(range(len(retention_pct)))
    ax2.set_xticklabels(retention_pct.index)
    ax2.grid(axis='y', alpha=0.3)
    
    # Add percentage values on bars
    for bar, pct, chrom in zip(bars_retention, retention_pct.values, retention_pct.index):
        height = bar.get_height()
        ax2.text(bar.get_x() + bar.get_width()/2., height + 1,
                f'{pct:.1f}%', ha='center', va='bottom', fontweight='bold')
        
        # Add raw and pruned counts below
        raw_count = raw_counts[chrom]
        pruned_count = pruned_counts[chrom]
        ax2.text(bar.get_x() + bar.get_width()/2., -5,
                f'Raw: {int(raw_count):,}\nPruned: {int(pruned_count):,}', 
                ha='center', va='top', fontsize=8)
    
    # Adjust ylim to make room for the text below bars
    ax2.set_ylim(bottom=-max(retention_pct.values)*0.15)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_SNP_Retention_By_Chromosome.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    plot_snp_distribution_main_chromosomes(data)
    plot_chromosome_comparison(data)
    plot_snp_retention_by_chromosome(data)