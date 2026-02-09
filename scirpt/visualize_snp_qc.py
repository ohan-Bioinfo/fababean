import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read the data
stages = ['01_Raw', '02_GENO', '02_MAF', '03_LD_Prune']
snps = [550755, 33077, 33077, 17170]
samples = [21, 21, 21, 21]
colors = ['#2E86AB', '#A23B72', '#F18F01', '#C73E1D']

# Create figure with subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: SNP Counts by QC Stage
bars = ax1.bar(stages, snps, color=colors, alpha=0.8, edgecolor='black')
ax1.set_title('SNP Counts by QC Stage\nFaba Bean GBS Data', fontsize=14, fontweight='bold')
ax1.set_ylabel('Number of SNPs', fontsize=12)
ax1.set_xlabel('QC Stage', fontsize=12)
ax1.tick_params(axis='x', rotation=45)

# Add value labels on bars
for bar, value in zip(bars, snps):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 1000,
             f'{value:,}', ha='center', va='bottom', fontweight='bold')

# Add retention percentages
retention_pct = [100, (33077/550755)*100, (33077/33077)*100, (17170/33077)*100]
for i, (stage, pct) in enumerate(zip(stages, retention_pct)):
    ax1.text(i, snps[i]//2, f'{pct:.1f}%', ha='center', va='center', 
             fontweight='bold', color='white', fontsize=10)

# Plot 2: SNP Retention Flow
ax2.plot(stages, snps, marker='o', linewidth=3, markersize=8, 
         color='#2E86AB', markerfacecolor='#A23B72')
ax2.set_title('SNP Retention Through QC Pipeline', fontsize=14, fontweight='bold')
ax2.set_ylabel('Number of SNPs (log scale)', fontsize=12)
ax2.set_xlabel('QC Stage', fontsize=12)
ax2.set_yscale('log')
ax2.tick_params(axis='x', rotation=45)
ax2.grid(True, alpha=0.3)

# Add value annotations
for i, (stage, count) in enumerate(zip(stages, snps)):
    ax2.annotate(f'{count:,}', (stage, count), textcoords="offset points", 
                 xytext=(0,10), ha='center', fontweight='bold')

plt.tight_layout()
plt.savefig('snp_qc_pipeline.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Pipeline visualization saved as 'snp_qc_pipeline.png'")
