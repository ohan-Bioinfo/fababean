import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import seaborn as sns

# Set style for clean visualization
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 11

def create_comprehensive_snp_heatmap():
    print("=== Creating Comprehensive SNP Density Heatmap ===")
    
    # Read the BIM file
    print("Reading SNP data from BIM file...")
    snp_data = pd.read_csv('Faba_chrOnly_pruned.bim', 
                          sep='\t', 
                          header=None, 
                          names=['chr', 'snp_id', 'cM', 'position', 'allele1', 'allele2'])
    
    print(f"Total SNPs loaded: {len(snp_data):,}")
    
    # Convert chromosome names and handle all classifications
    snp_data['chr'] = snp_data['chr'].astype(str)
    
    # Standardize chromosome names - add 'chr' prefix to numeric ones
    def standardize_chrom_name(chrom):
        if chrom.isdigit():
            return f'chr{chrom}'
        elif chrom in ['chr1L', 'chr1S']:
            return chrom
        else:
            return f'chr{chrom}'  # Handle other cases
    
    snp_data['chr_standardized'] = snp_data['chr'].apply(standardize_chrom_name)
    
    # Get all unique chromosomes after standardization
    all_chroms = sorted(snp_data['chr_standardized'].unique())
    print(f"All chromosome classifications: {all_chroms}")
    print(f"Total chromosome types: {len(all_chroms)}")
    
    # Create density matrix
    bin_size = 1000000  # 1Mb
    max_bins = 0
    chrom_lengths = {}
    
    print("\nCalculating chromosome lengths and SNP distributions...")
    for chrom in all_chroms:
        chrom_data = snp_data[snp_data['chr_standardized'] == chrom]
        if len(chrom_data) > 0:
            max_pos = chrom_data['position'].max()
            chrom_lengths[chrom] = max_pos
            num_bins = int(np.ceil(max_pos / bin_size))
            max_bins = max(max_bins, num_bins)
            print(f"  {chrom}: {len(chrom_data):,} SNPs, {num_bins} bins, max pos: {max_pos:,} bp")

    # Create the density matrix (chromosomes as rows, bins as columns)
    density_matrix = np.zeros((len(all_chroms), max_bins))
    
    for i, chrom in enumerate(all_chroms):
        chrom_data = snp_data[snp_data['chr_standardized'] == chrom]
        if len(chrom_data) > 0:
            max_pos = chrom_lengths[chrom]
            bins = np.arange(0, max_pos + bin_size, bin_size)
            counts, _ = np.histogram(chrom_data['position'], bins=bins)
            
            # Fill the row for this chromosome
            for j in range(min(len(counts), max_bins)):
                density_matrix[i, j] = counts[j]
    
    # Create the heatmap figure
    fig, ax = plt.subplots(figsize=(16, 10))
    
    # Create heatmap with better color scaling
    vmax = np.percentile(density_matrix[density_matrix > 0], 95) if np.any(density_matrix > 0) else 1
    heatmap = sns.heatmap(density_matrix, 
                         cmap='YlOrRd',
                         cbar_kws={'label': 'SNP Count per 1Mb Bin', 'shrink': 0.8},
                         ax=ax,
                         square=False,
                         vmin=0,
                         vmax=vmax)
    
    # Set y-axis labels (all chromosomes on the left)
    ax.set_yticks(np.arange(len(all_chroms)) + 0.5)
    ax.set_yticklabels(all_chroms, rotation=0, fontsize=11, fontweight='bold')
    ax.set_ylabel('Chromosome', fontsize=13, fontweight='bold', labelpad=15)
    
    # Set x-axis labels (bins) - show fewer labels to avoid overlap
    if max_bins > 50:
        x_tick_interval = max(1, max_bins // 15)
    else:
        x_tick_interval = max(1, max_bins // 10)
    
    x_ticks = np.arange(0, max_bins, x_tick_interval)
    x_labels = [f'{x}' for x in x_ticks]
    ax.set_xticks(x_ticks)
    ax.set_xticklabels(x_labels, rotation=0, fontsize=9)
    ax.set_xlabel('Genomic Position (1Mb Bins)', fontsize=13, fontweight='bold', labelpad=10)
    
    # Set title
    total_snps = len(snp_data)
    ax.set_title(f'Comprehensive SNP Density Distribution\n{total_snps:,} SNPs across {len(all_chroms)} chromosomes | 21 Faba Bean Accessions', 
                 fontsize=16, fontweight='bold', pad=20)
    
    # Improve layout
    plt.tight_layout()
    
    # Save high-quality output
    plt.savefig('SNP_density_comprehensive_heatmap.png', dpi=300, bbox_inches='tight')
    plt.savefig('SNP_density_comprehensive_heatmap.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Print comprehensive summary
    print(f"\n📊 COMPREHENSIVE SUMMARY:")
    print(f"Total SNPs analyzed: {total_snps:,}")
    print(f"Chromosome types: {len(all_chroms)}")
    print(f"Matrix dimensions: {density_matrix.shape[0]} chromosomes × {density_matrix.shape[1]} bins")
    print(f"Bin size: 1 Mb")
    print(f"Accessions: 21")
    
    # Detailed chromosome statistics
    print(f"\n📈 DETAILED CHROMOSOME STATISTICS:")
    for i, chrom in enumerate(all_chroms):
        chrom_snps = np.sum(density_matrix[i, :])
        if chrom_snps > 0:
            max_density = np.max(density_matrix[i, :])
            non_zero_bins = np.sum(density_matrix[i, :] > 0)
            mean_density = np.mean(density_matrix[i, density_matrix[i, :] > 0])
            coverage = (non_zero_bins / max_bins) * 100
            print(f"  {chrom:8}: {chrom_snps:5.0f} SNPs, max: {max_density:3.0f}/Mb, mean: {mean_density:4.1f}/Mb, coverage: {coverage:4.1f}%")
    
    # Overall statistics
    total_non_zero_bins = np.sum(density_matrix > 0)
    total_possible_bins = density_matrix.shape[0] * density_matrix.shape[1]
    overall_coverage = (total_non_zero_bins / total_possible_bins) * 100
    
    print(f"\n🌍 OVERALL STATISTICS:")
    print(f"Total genomic bins: {total_possible_bins:,}")
    print(f"Bins with SNPs: {total_non_zero_bins:,} ({overall_coverage:.1f}% coverage)")
    print(f"Average SNP density: {total_snps/(np.sum(list(chrom_lengths.values()))/1e6):.1f} SNPs/Mb")
    
    # Save detailed data
    density_df = pd.DataFrame(density_matrix, 
                             index=all_chroms,
                             columns=[f'Bin_{i}' for i in range(max_bins)])
    density_df.to_csv('SNP_density_comprehensive_matrix.csv')
    print(f"\n💾 Detailed matrix saved as 'SNP_density_comprehensive_matrix.csv'")

def create_filtered_heatmap():
    """Create heatmap focusing on main chromosomes only"""
    print("\n=== Creating Filtered Heatmap (Main Chromosomes Only) ===")
    
    # Read and process data
    snp_data = pd.read_csv('Faba_chrOnly_pruned.bim', 
                          sep='\t', 
                          header=None, 
                          names=['chr', 'snp_id', 'cM', 'position', 'allele1', 'allele2'])
    
    snp_data['chr'] = snp_data['chr'].astype(str)
    
    # Focus on main chromosomes only
    main_chroms = ['chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']
    # Add 'chr' prefix to any numeric chromosomes in the main list
    main_chroms_extended = main_chroms + [f'chr{x}' for x in range(1, 7) if f'chr{x}' not in main_chroms]
    
    # Filter for main chromosomes
    snp_data_main = snp_data[snp_data['chr'].isin(main_chroms_extended) | 
                            snp_data['chr'].apply(lambda x: f'chr{x}' if x.isdigit() else x).isin(main_chroms_extended)]
    
    print(f"SNPs on main chromosomes: {len(snp_data_main):,}")
    
    # Standardize names for main chromosomes
    def get_main_chrom_name(chrom):
        chrom_str = str(chrom)
        if chrom_str.isdigit():
            return f'chr{chrom_str}'
        elif chrom_str in ['1L', '1S']:
            return f'chr{chrom_str}'
        else:
            return chrom_str
    
    snp_data_main['chrom_clean'] = snp_data_main['chr'].apply(get_main_chrom_name)
    main_chroms_clean = sorted(snp_data_main['chrom_clean'].unique())
    
    print(f"Main chromosomes found: {main_chroms_clean}")
    
    # Create density matrix for main chromosomes
    bin_size = 1000000
    max_bins = 0
    chrom_lengths_main = {}
    
    for chrom in main_chroms_clean:
        chrom_data = snp_data_main[snp_data_main['chrom_clean'] == chrom]
        if len(chrom_data) > 0:
            max_pos = chrom_data['position'].max()
            chrom_lengths_main[chrom] = max_pos
            num_bins = int(np.ceil(max_pos / bin_size))
            max_bins = max(max_bins, num_bins)
    
    density_matrix_main = np.zeros((len(main_chroms_clean), max_bins))
    
    for i, chrom in enumerate(main_chroms_clean):
        chrom_data = snp_data_main[snp_data_main['chrom_clean'] == chrom]
        if len(chrom_data) > 0:
            max_pos = chrom_lengths_main[chrom]
            bins = np.arange(0, max_pos + bin_size, bin_size)
            counts, _ = np.histogram(chrom_data['position'], bins=bins)
            density_matrix_main[i, :len(counts)] = counts
    
    # Create the filtered heatmap
    plt.figure(figsize=(14, 8))
    
    vmax = np.percentile(density_matrix_main[density_matrix_main > 0], 95) if np.any(density_matrix_main > 0) else 1
    plt.imshow(density_matrix_main, cmap='YlOrRd', aspect='auto', vmin=0, vmax=vmax)
    
    plt.yticks(range(len(main_chroms_clean)), main_chroms_clean, fontsize=12, fontweight='bold')
    
    # X-axis with reasonable label spacing
    x_ticks = range(0, max_bins, max(1, max_bins//10))
    plt.xticks(x_ticks, [str(x) for x in x_ticks], fontsize=10)
    
    plt.xlabel('Genomic Position (1Mb Bins)', fontsize=13, fontweight='bold', labelpad=10)
    plt.ylabel('Chromosome', fontsize=13, fontweight='bold', labelpad=15)
    
    cbar = plt.colorbar(shrink=0.8)
    cbar.set_label('SNP Count per 1Mb Bin', rotation=270, labelpad=20)
    
    plt.title(f'SNP Density: Main Chromosomes\n{len(snp_data_main):,} SNPs, 21 Accessions', 
              fontsize=16, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig('SNP_density_main_chromosomes.png', dpi=300, bbox_inches='tight')
    plt.show()
    
    print(f"✓ Filtered heatmap created with {len(snp_data_main):,} SNPs on {len(main_chroms_clean)} main chromosomes")

if __name__ == "__main__":
    # Run comprehensive analysis (all chromosome classifications)
    create_comprehensive_snp_heatmap()
    
    # Also create filtered version focusing on main chromosomes
    create_filtered_heatmap()