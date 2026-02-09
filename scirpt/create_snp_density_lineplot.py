import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read and clean the data
bim_file = "03_LD_Prune/Faba_chrOnly_pruned.bim"
df = pd.read_csv(bim_file, sep='\t', header=None, 
                 names=['chr', 'snp_id', 'cM', 'pos', 'a1', 'a2'])

# Fix chromosome names
df['chr'] = df['chr'].apply(lambda x: f"chr{x}" if str(x).isdigit() else x)
main_chromosomes = ['chr1L', 'chr1S', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6']
df = df[df['chr'].isin(main_chromosomes)]

# Remove position outliers
df = df[df['pos'] <= 500000000]

# Create the plot
fig, ax = plt.subplots(figsize=(14, 8))

colors = plt.cm.Set3(np.linspace(0, 1, len(main_chromosomes)))

for i, chrom in enumerate(main_chromosomes):
    chrom_data = df[df['chr'] == chrom]
    if len(chrom_data) > 0:
        # Create density using rolling window
        positions = chrom_data['pos'].sort_values().values
        if len(positions) > 1:
            # Simple density: count SNPs in 1Mb windows
            window_size = 1000000
            max_pos = positions.max()
            windows = np.arange(0, max_pos + window_size, window_size)
            densities = []
            midpoints = []
            
            for j in range(len(windows)-1):
                count = np.sum((positions >= windows[j]) & (positions < windows[j+1]))
                densities.append(count)
                midpoints.append((windows[j] + windows[j+1]) / 2)
            
            # Plot as line
            ax.plot(midpoints, densities, label=chrom, color=colors[i], linewidth=2)
            ax.fill_between(midpoints, 0, densities, alpha=0.3, color=colors[i])

ax.set_title('SNP Density Across Chromosomes\n1Mb Rolling Windows', 
             fontsize=16, fontweight='bold')
ax.set_xlabel('Genomic Position (bp)', fontsize=12)
ax.set_ylabel('SNP Count per 1Mb Window', fontsize=12)
ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('Faba_Bean_SNP_Density_Lineplot.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ SNP density line plot saved as 'Faba_Bean_SNP_Density_Lineplot.png'")
