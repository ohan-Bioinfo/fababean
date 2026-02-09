import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read the BIM file
bim_file = "03_LD_Prune/Faba_chrOnly_pruned.bim"
df = pd.read_csv(bim_file, sep='\t', header=None, 
                 names=['chr', 'snp_id', 'cM', 'pos', 'a1', 'a2'])

# Filter only main chromosomes
main_chromosomes = ['chr1L', 'chr1S', '2', '3', '4', '5', '6']
df = df[df['chr'].isin(main_chromosomes)]

# Define chromosome order for plotting
chr_order = ['chr1L', 'chr1S', '2', '3', '4', '5', '6']

# Create 5Mb windows and count SNPs
window_size = 5000000  # 5Mb in base pairs

# Get chromosome lengths (approximate from max positions)
chr_lengths = {}
for chrom in chr_order:
    chrom_data = df[df['chr'] == chrom]
    if len(chrom_data) > 0:
        chr_lengths[chrom] = chrom_data['pos'].max()
    else:
        chr_lengths[chrom] = 0

print("Chromosome lengths (approximate):")
for chrom, length in chr_lengths.items():
    print(f"  {chrom}: {length:,} bp")

# Create windows and count SNPs
heatmap_data = []
for chrom in chr_order:
    max_pos = chr_lengths[chrom]
    if max_pos == 0:
        continue
        
    # Create windows
    windows = list(range(0, max_pos + window_size, window_size))
    for i in range(len(windows) - 1):
        start = windows[i]
        end = windows[i + 1]
        
        # Count SNPs in this window
        snp_count = len(df[(df['chr'] == chrom) & 
                          (df['pos'] >= start) & 
                          (df['pos'] < end)])
        
        heatmap_data.append({
            'chromosome': chrom,
            'window_start': start,
            'window_end': end,
            'snp_count': snp_count,
            'window_mid': (start + end) / 2,
            'window_id': f"{chrom}_{i}"
        })

# Convert to DataFrame
heatmap_df = pd.DataFrame(heatmap_data)

# Pivot for heatmap
pivot_df = heatmap_df.pivot_table(
    index='chromosome', 
    columns='window_mid', 
    values='snp_count',
    fill_value=0
)

# Reindex to maintain chromosome order
pivot_df = pivot_df.reindex(chr_order)

# Create the heatmap
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 10))

# Plot 1: Heatmap
im = ax1.imshow(pivot_df.values, aspect='auto', cmap='YlOrRd', 
                interpolation='nearest')
ax1.set_yticks(range(len(pivot_df.index)))
ax1.set_yticklabels(pivot_df.index)
ax1.set_title('SNP Density Heatmap - 5Mb Windows\nFaba Bean GBS (LD-pruned SNPs)', 
              fontsize=16, fontweight='bold')
ax1.set_ylabel('Chromosome', fontsize=12)

# Add colorbar
cbar = plt.colorbar(im, ax=ax1)
cbar.set_label('SNP Count per 5Mb Window', fontsize=12)

# Add values to heatmap cells
for i in range(len(pivot_df.index)):
    for j in range(len(pivot_df.columns)):
        value = pivot_df.iloc[i, j]
        if value > 0:  # Only label non-zero cells
            ax1.text(j, i, f'{int(value)}', 
                    ha='center', va='center', 
                    fontsize=8, fontweight='bold',
                    color='black' if value < np.max(pivot_df.values) * 0.7 else 'white')

# Plot 2: Bar plot of total SNPs per chromosome
chrom_totals = df['chr'].value_counts().reindex(chr_order)
bars = ax2.bar(range(len(chrom_totals)), chrom_totals.values, 
               color=plt.cm.Set3(np.linspace(0, 1, len(chrom_totals))))
ax2.set_title('Total SNPs per Chromosome', fontsize=14, fontweight='bold')
ax2.set_xlabel('Chromosome', fontsize=12)
ax2.set_ylabel('Number of SNPs', fontsize=12)
ax2.set_xticks(range(len(chrom_totals)))
ax2.set_xticklabels(chrom_totals.index)

# Add values on bars
for bar, value in zip(bars, chrom_totals.values):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 10,
            f'{value:,}', ha='center', va='bottom', 
            fontweight='bold', fontsize=10)

plt.tight_layout()
plt.savefig('Faba_Bean_SNP_Density_Heatmap_5Mb.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ SNP density heatmap saved as 'Faba_Bean_SNP_Density_Heatmap_5Mb.png'")

# Save the heatmap data
heatmap_df.to_csv('snp_density_5mb_windows.csv', index=False)
print("✓ SNP density data saved as 'snp_density_5mb_windows.csv'")

# Print summary statistics
print(f"\n📊 SNP DENSITY SUMMARY:")
print(f"Total windows analyzed: {len(heatmap_df)}")
print(f"Average SNPs per 5Mb window: {heatmap_df['snp_count'].mean():.1f}")
print(f"Maximum SNPs in a 5Mb window: {heatmap_df['snp_count'].max()}")
print(f"Windows with zero SNPs: {len(heatmap_df[heatmap_df['snp_count'] == 0])}")

# Print per-chromosome statistics
print(f"\n📈 PER-CHROMOSOME DENSITY:")
for chrom in chr_order:
    chrom_data = heatmap_df[heatmap_df['chromosome'] == chrom]
    if len(chrom_data) > 0:
        avg_density = chrom_data['snp_count'].mean()
        max_density = chrom_data['snp_count'].max()
        print(f"  {chrom}: {avg_density:.1f} SNPs/5Mb (max: {max_density})")
