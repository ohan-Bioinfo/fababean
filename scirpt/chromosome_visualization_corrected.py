import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Create figure
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

# Data for pipeline stages
stages = ['01_Raw', '02_GENO', '02_MAF', '03_LD_Prune']
snps = [550755, 33077, 33077, 17170]
retention_pct = [100, 6.0, 100, 51.9]
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

# Plot 1: SNP QC Pipeline
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

# Plot 2: Chromosome Distribution with your specific names and order
import subprocess

# Get chromosome counts
result = subprocess.run(["cut", "-f1", "03_LD_Prune/Faba_chrOnly_pruned.bim"], 
                       capture_output=True, text=True)
chromosomes = result.stdout.strip().split('\n')

# Count SNPs per chromosome
chrom_counts = {}
for chrom in chromosomes:
    chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

# Define the specific order you want
desired_order = ['chr1L', 'chr1S', '2', '3', '4', '5', '6']

# Filter and sort according to desired order
sorted_chroms = []
for chrom in desired_order:
    if chrom in chrom_counts:
        sorted_chroms.append((chrom, chrom_counts[chrom]))

chrom_names = [x[0] for x in sorted_chroms]
chrom_snps = [x[1] for x in sorted_chroms]

# Create colors for each chromosome
colors_chrom = plt.cm.Set3(np.linspace(0, 1, len(chrom_names)))
bars_chrom = ax2.bar(chrom_names, chrom_snps, color=colors_chrom, alpha=0.8, edgecolor='black')
ax2.set_title('SNP Distribution Across Chromosomes\n(Final LD-pruned Dataset)', 
             fontsize=14, fontweight='bold')
ax2.set_ylabel('Number of SNPs', fontsize=12)
ax2.set_xlabel('Chromosome', fontsize=12)

# Add values on top of bars
total_snps = sum(chrom_snps)
for bar, value, chrom in zip(bars_chrom, chrom_snps, chrom_names):
    height = bar.get_height()
    percentage = (value / total_snps) * 100
    ax2.text(bar.get_x() + bar.get_width()/2., height + 20,
            f'{value:,}', ha='center', va='bottom', 
            fontweight='bold', fontsize=11)
    # Add percentage inside the bar
    ax2.text(bar.get_x() + bar.get_width()/2., height/2,
            f'{percentage:.1f}%', ha='center', va='center', 
            fontweight='bold', color='white', fontsize=10)

plt.tight_layout()
plt.savefig('Faba_Bean_Chromosomes_Corrected.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Corrected chromosome visualization saved as 'Faba_Bean_Chromosomes_Corrected.png'")

# Print summary
print(f"\n📊 CHROMOSOME SNP DISTRIBUTION (Your Order):")
for chrom, count in sorted_chroms:
    pct = (count / total_snps) * 100
    print(f"  {chrom}: {count:,} SNPs ({pct:.1f}%)")
print(f"  Total SNPs: {total_snps:,}")

# Create separate heterozygosity plot
fig, ax = plt.subplots(figsize=(12, 6))
het_data = pd.read_csv('heterozygosity_summary.csv')

bars_het = ax.bar(het_data['IID'], het_data['het_rate'], 
                  color='#4ECDC4', alpha=0.8, edgecolor='black')
ax.axhline(y=het_data['het_rate'].mean(), color='red', linestyle='--', 
           linewidth=2, label=f'Mean: {het_data["het_rate"].mean():.3f}')
ax.set_title('Sample Heterozygosity Rates - No Outliers Detected', 
            fontsize=14, fontweight='bold')
ax.set_ylabel('Heterozygosity Rate', fontsize=12)
ax.set_xlabel('Sample ID', fontsize=12)
ax.legend()
ax.tick_params(axis='x', rotation=45)

# Add heterozygosity values on bars
for bar, rate in zip(bars_het, het_data['het_rate']):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 0.005,
            f'{rate:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('Faba_Bean_Heterozygosity.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Heterozygosity visualization saved as 'Faba_Bean_Heterozygosity.png'")
