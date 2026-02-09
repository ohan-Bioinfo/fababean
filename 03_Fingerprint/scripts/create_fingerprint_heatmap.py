# scripts/create_fingerprint_heatmap.py

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

def load_genotype_data():
    """Load genotype data from PLINK files"""
    # Read sample information
    fam = pd.read_csv('data/faba_fingerprint.fam',
                      delim_whitespace=True, header=None,
                      names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
    
    # Read SNP information
    bim = pd.read_csv('data/faba_fingerprint.bim',
                      delim_whitespace=True, header=None,
                      names=['CHR', 'SNP', 'cM', 'POS', 'A1', 'A2'])
    
    # Read genotype data (from .ped file)
    ped = pd.read_csv('output/faba_fingerprint_genotypes.ped',
                      delim_whitespace=True, header=None)
    
    # Extract genotype calls (columns 6 onwards are genotypes, pairs of alleles)
    genotype_cols = ped.iloc[:, 6:]
    sample_ids = fam['IID'].astype(str).tolist()

    # Convert genotype pairs to single genotype calls
    genotypes = []
    for i in range(0, genotype_cols.shape[1], 2):
        snp_genotypes = []
        for j in range(len(genotype_cols)):
            allele1 = genotype_cols.iloc[j, i]
            allele2 = genotype_cols.iloc[j, i + 1]
            if allele1 == '0' or allele2 == '0':
                snp_genotypes.append('Missing')
            elif allele1 == allele2:
                snp_genotypes.append(f'{allele1}/{allele1}')
            else:
                snp_genotypes.append(f'{allele1}/{allele2}')
        genotypes.append(snp_genotypes)
    
    # Create genotype matrix with SNP1, SNP2, ..., SNP150
    snp_labels = [f'SNP{i+1}' for i in range(len(genotypes))]
    genotype_df = pd.DataFrame(genotypes, index=snp_labels, columns=sample_ids).T

    # Sort accession/sample IDs numerically (as strings)
    def numeric_key(x):
        try:
            return int(x)
        except ValueError:
            return float('inf')
    
    genotype_df = genotype_df.sort_index(key=lambda x: [numeric_key(i) for i in x])

    return genotype_df, bim, fam

def create_heatmap(genotype_df, output_file):
    """Create categorical heatmap of genotypes"""
    # Define color mapping
    color_map = {
        'A/A': 'blue',
        'C/C': 'green',
        'G/G': 'red',
        'T/T': 'purple',
        'Missing': 'white'
    }

    # Define heterozygous genotypes
    heterozygous_genotypes = [
        'A/C', 'A/G', 'A/T', 'C/G', 'C/T', 'G/T',
        'C/A', 'G/A', 'T/A', 'G/C', 'T/C', 'T/G'
    ]
    for het in heterozygous_genotypes:
        color_map[het] = 'gray'
    
    # Convert to color matrix
    color_matrix = genotype_df.applymap(lambda x: color_map.get(x, 'white'))

    num_accessions, num_snps = color_matrix.shape

    # Adjust figure size based on SNP count for compact layout
    fig_width = max(10, num_snps * 0.15)
    fig_height = max(8, num_accessions * 0.2)

    fig, ax = plt.subplots(figsize=(fig_width, fig_height))

    # Draw each cell
    for i in range(num_accessions):
        for j in range(num_snps):
            ax.add_patch(plt.Rectangle(
                (j, i), 1, 1,
                facecolor=color_matrix.iloc[i, j],
                edgecolor='white', lw=0.3
            ))

    # Set ticks
    ax.set_xticks(np.arange(num_snps) + 0.5)
    ax.set_yticks(np.arange(num_accessions) + 0.5)

    # X-axis: SNP labels (number bold only), smaller font
    xtick_labels = [rf'SNP$\bf{{{i+1}}}$' for i in range(num_snps)]
    ax.set_xticklabels(xtick_labels, rotation=90, fontsize=5)

    # Y-axis: sorted accession IDs
    ax.set_yticklabels(genotype_df.index, fontsize=7, fontweight='bold')

    ax.set_xlim(0, num_snps)
    ax.set_ylim(0, num_accessions)
    ax.set_aspect('equal')

    # Legend
    legend_elements = [
        Patch(facecolor='blue', label='A/A'),
        Patch(facecolor='green', label='C/C'),
        Patch(facecolor='red', label='G/G'),
        Patch(facecolor='purple', label='T/T'),
        Patch(facecolor='gray', label='Heterozygotes'),
        Patch(facecolor='white', label='Missing')
    ]
    ax.legend(handles=legend_elements, loc='upper right',
              bbox_to_anchor=(1.15, 1), fontsize=8)

    # Titles
    plt.title('Faba Bean Fingerprint Panel - 150 SNPs',
              fontsize=14, fontweight='bold', pad=20)
    plt.xlabel('SNP Markers', fontsize=10, fontweight='bold')
    plt.ylabel('Accessions', fontsize=10, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.savefig(output_file.replace('.png', '.pdf'), bbox_inches='tight')
    plt.show()

    return color_map

def main():
    print("Loading genotype data...")
    genotype_df, bim, fam = load_genotype_data()

    print(f"Genotype matrix shape: {genotype_df.shape}")
    print(f"Samples: {len(genotype_df)}")
    print(f"SNPs: {len(genotype_df.columns)}")

    print("Creating fingerprint heatmap...")
    color_map = create_heatmap(genotype_df, 'plots/faba_fingerprint_heatmap_categorical.png')

    # Save summary to CSV
    genotype_summary = genotype_df.apply(pd.Series.value_counts).fillna(0)
    genotype_summary.to_csv('output/genotype_summary.csv')

    print("Fingerprint analysis complete!")
    print("Output files saved in plots/ and output/ directories")

if __name__ == "__main__":
    main()
