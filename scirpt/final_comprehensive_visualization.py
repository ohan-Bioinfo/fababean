import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Create figure with multiple subplots
fig = plt.figure(figsize=(16, 12))

# Define grid layout
gs = plt.GridSpec(3, 2, figure=fig)

# Data
stages = ['01_Raw', '02_GENO', '02_MAF', '03_LD_Prune']
snps = [550755, 33077, 33077, 17170]
retention_pct = [100, 6.0, 100, 51.9]
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

# Plot 1: SNP QC Pipeline (top left)
ax1 = fig.add_subplot(gs[0, 0])
bars = ax1.bar(stages, snps, color=colors, alpha=0.8, edgecolor='black')
ax1.set_title('SNP QC Pipeline - Faba Bean GBS', fontsize=14, fontweight='bold')
ax1.set_ylabel('Number of SNPs', fontsize=12)
ax1.set_yscale('log')

# Add value labels
for bar, value, pct in zip(bars, snps, retention_pct):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height * 1.1,
             f'{value:,}\n({pct}%)', ha='center', va='bottom', 
             fontweight='bold', fontsize=9)

# Plot 2: Retention Flow (top right)
ax2 = fig.add_subplot(gs[0, 1])
ax2.plot(stages, snps, marker='o', linewidth=3, markersize=10, 
         color='#2E86AB', markerfacecolor='#A23B72')
ax2.set_title('SNP Retention Flow', fontsize=14, fontweight='bold')
ax2.set_ylabel('SNP Count (log scale)', fontsize=12)
ax2.set_yscale('log')
ax2.grid(True, alpha=0.3)

# Add retention percentages
for i, (stage, count, pct) in enumerate(zip(stages, snps, retention_pct)):
    if i > 0:
        ax2.annotate(f'{pct}%', (stage, count), textcoords="offset points", 
                    xytext=(0,15), ha='center', fontweight='bold', fontsize=10)

# Plot 3: Chromosome Distribution (middle left)
# Get chromosome counts
import subprocess
result = subprocess.run(["cut", "-f1", "03_LD_Prune/Faba_chrOnly_pruned.bim"], 
                       capture_output=True, text=True)
chromosomes = result.stdout.strip().split('\n')
chrom_counts = {}
for chrom in chromosomes:
    chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

# Sort chromosomes
sorted_chroms = sorted(chrom_counts.items(), key=lambda x: int(x[0][3:]) if x[0].startswith('chr') else x[0])
chrom_names = [x[0] for x in sorted_chroms]
chrom_snps = [x[1] for x in sorted_chroms]

ax3 = fig.add_subplot(gs[1, :])
colors_chrom = plt.cm.Set3(np.linspace(0, 1, len(chrom_names)))
bars_chrom = ax3.bar(chrom_names, chrom_snps, color=colors_chrom, alpha=0.8, edgecolor='black')
ax3.set_title('SNP Distribution Across Chromosomes (Final LD-pruned Dataset)', 
             fontsize=14, fontweight='bold')
ax3.set_ylabel('Number of SNPs', fontsize=12)
ax3.set_xlabel('Chromosome', fontsize=12)

# Add values on bars
total_snps = sum(chrom_snps)
for bar, value in zip(bars_chrom, chrom_snps):
    height = bar.get_height()
    percentage = (value / total_snps) * 100
    ax3.text(bar.get_x() + bar.get_width()/2., height + 10,
            f'{value:,}\n({percentage:.1f}%)', ha='center', va='bottom', 
            fontweight='bold', fontsize=9)

# Plot 4: Heterozygosity Distribution (bottom)
het_data = pd.read_csv('heterozygosity_summary.csv')
ax4 = fig.add_subplot(gs[2, :])
sample_order = het_data.sort_values('het_rate')['IID']
colors_het = ['#4ECDC4' for _ in range(len(het_data))]  # All green - no outliers

bars_het = ax4.bar(het_data['IID'], het_data['het_rate'], 
                   color=colors_het, alpha=0.8, edgecolor='black')
ax4.axhline(y=het_data['het_rate'].mean(), color='red', linestyle='--', 
            linewidth=2, label=f'Mean: {het_data["het_rate"].mean():.3f}')
ax4.set_title('Sample Heterozygosity Rates - No Outliers Detected', 
             fontsize=14, fontweight='bold')
ax4.set_ylabel('Heterozygosity Rate', fontsize=12)
ax4.set_xlabel('Sample ID', fontsize=12)
ax4.legend()
ax4.tick_params(axis='x', rotation=45)

# Add heterozygosity values on bars
for bar, rate in zip(bars_het, het_data['het_rate']):
    height = bar.get_height()
    ax4.text(bar.get_x() + bar.get_width()/2., height + 0.005,
             f'{rate:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('Faba_Bean_Complete_QC_Summary.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Complete QC summary visualization saved as 'Faba_Bean_Complete_QC_Summary.png'")
