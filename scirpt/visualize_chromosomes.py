import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Get chromosome counts from the pruned BIM file
import subprocess
result = subprocess.run(["cut", "-f1", "03_LD_Prune/Faba_chrOnly_pruned.bim"], 
                       capture_output=True, text=True)
chromosomes = result.stdout.strip().split('\n')

# Count SNPs per chromosome
chrom_counts = {}
for chrom in chromosomes:
    chrom_counts[chrom] = chrom_counts.get(chrom, 0) + 1

# Sort by chromosome number
sorted_chroms = sorted(chrom_counts.items(), key=lambda x: int(x[0][3:]) if x[0].startswith('chr') else x[0])
chrom_names = [x[0] for x in sorted_chroms]
chrom_snps = [x[1] for x in sorted_chroms]

# Create color palette
colors = plt.cm.Set3(np.linspace(0, 1, len(chrom_names)))

# Create plot
fig, ax = plt.subplots(figsize=(12, 6))
bars = ax.bar(chrom_names, chrom_snps, color=colors, alpha=0.8, edgecolor='black')

# Customize plot
ax.set_title('SNP Distribution Across Chromosomes\nFaba Bean - LD Pruned Dataset', 
             fontsize=14, fontweight='bold')
ax.set_ylabel('Number of SNPs', fontsize=12)
ax.set_xlabel('Chromosome', fontsize=12)

# Add value labels on top of bars
for bar, value in zip(bars, chrom_snps):
    height = bar.get_height()
    ax.text(bar.get_x() + bar.get_width()/2., height + 20,
            f'{value:,}', ha='center', va='bottom', fontweight='bold', fontsize=10)

# Add percentage labels inside bars
total_snps = sum(chrom_snps)
for bar, value in zip(bars, chrom_snps):
    height = bar.get_height()
    percentage = (value / total_snps) * 100
    ax.text(bar.get_x() + bar.get_width()/2., height/2,
            f'{percentage:.1f}%', ha='center', va='center', 
            fontweight='bold', color='black', fontsize=9)

plt.xticks(rotation=45)
plt.tight_layout()
plt.savefig('chromosome_snp_distribution.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Chromosome distribution saved as 'chromosome_snp_distribution.png'")

# Print summary
print(f"\n📊 CHROMOSOME SNP SUMMARY:")
print(f"Total SNPs: {total_snps:,}")
for chrom, count in sorted_chroms:
    pct = (count / total_snps) * 100
    print(f"{chrom}: {count:,} SNPs ({pct:.1f}%)")
