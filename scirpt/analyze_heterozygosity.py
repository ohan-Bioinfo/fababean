import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib import rcParams
import seaborn as sns

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read heterozygosity data
het_data = pd.read_csv('04_Het_QC/Faba_chrOnly_het.het', sep='\s+')

# Calculate heterozygosity rate
het_data['het_rate'] = (het_data['N(NM)'] - het_data['O(HOM)']) / het_data['N(NM)']

# Calculate statistics
mean_het = het_data['het_rate'].mean()
std_het = het_data['het_rate'].std()
threshold_high = mean_het + (3 * std_het)
threshold_low = mean_het - (3 * std_het)

# Identify outliers
outliers_high = het_data[het_data['het_rate'] > threshold_high]
outliers_low = het_data[het_data['het_rate'] < threshold_low]

print("=== HETEROZYGOSITY ANALYSIS ===")
print(f"Mean heterozygosity: {mean_het:.4f}")
print(f"Standard deviation: {std_het:.4f}")
print(f"Upper threshold (mean + 3SD): {threshold_high:.4f}")
print(f"Lower threshold (mean - 3SD): {threshold_low:.4f}")
print(f"High heterozygosity outliers: {len(outliers_high)}")
print(f"Low heterozygosity outliers: {len(outliers_low)}")

if len(outliers_high) > 0:
    print("\nHigh heterozygosity outliers (possible contamination):")
    for _, row in outliers_high.iterrows():
        print(f"  {row['FID']} {row['IID']}: {row['het_rate']:.4f}")

if len(outliers_low) > 0:
    print("\nLow heterozygosity outliers (possible inbreeding/relatedness):")
    for _, row in outliers_low.iterrows():
        print(f"  {row['FID']} {row['IID']}: {row['het_rate']:.4f}")

# Create visualization
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Heterozygosity distribution with outliers
bars = ax1.bar(het_data['IID'], het_data['het_rate'], 
               color=['#FF6B6B' if x > threshold_high else '#4ECDC4' if x < threshold_low else '#45B7D1' for x in het_data['het_rate']],
               alpha=0.8, edgecolor='black')
ax1.axhline(y=mean_het, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_het:.4f}')
ax1.axhline(y=threshold_high, color='orange', linestyle='--', linewidth=1, label=f'+3SD: {threshold_high:.4f}')
ax1.axhline(y=threshold_low, color='orange', linestyle='--', linewidth=1, label=f'-3SD: {threshold_low:.4f}')
ax1.set_title('Sample Heterozygosity Rates\nFaba Bean GBS Data', fontsize=14, fontweight='bold')
ax1.set_ylabel('Heterozygosity Rate', fontsize=12)
ax1.set_xlabel('Sample ID', fontsize=12)
ax1.tick_params(axis='x', rotation=45)
ax1.legend()

# Add value labels on bars
for bar, rate in zip(bars, het_data['het_rate']):
    height = bar.get_height()
    ax1.text(bar.get_x() + bar.get_width()/2., height + 0.001,
             f'{rate:.3f}', ha='center', va='bottom', fontsize=8, fontweight='bold')

# Plot 2: Histogram of heterozygosity rates
ax2.hist(het_data['het_rate'], bins=10, color='#45B7D1', alpha=0.7, edgecolor='black')
ax2.axvline(mean_het, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_het:.4f}')
ax2.axvline(threshold_high, color='orange', linestyle='--', linewidth=1, label=f'+3SD: {threshold_high:.4f}')
ax2.axvline(threshold_low, color='orange', linestyle='--', linewidth=1, label=f'-3SD: {threshold_low:.4f}')
ax2.set_title('Distribution of Heterozygosity Rates', fontsize=14, fontweight='bold')
ax2.set_xlabel('Heterozygosity Rate', fontsize=12)
ax2.set_ylabel('Number of Samples', fontsize=12)
ax2.legend()

# Add outlier annotations
if len(outliers_high) > 0:
    for _, outlier in outliers_high.iterrows():
        ax2.annotate(f"High: {outlier['IID']}", 
                    xy=(outlier['het_rate'], 0), 
                    xytext=(outlier['het_rate'], 1),
                    arrowprops=dict(arrowstyle='->', color='red'),
                    fontweight='bold')

if len(outliers_low) > 0:
    for _, outlier in outliers_low.iterrows():
        ax2.annotate(f"Low: {outlier['IID']}", 
                    xy=(outlier['het_rate'], 0), 
                    xytext=(outlier['het_rate'], 1),
                    arrowprops=dict(arrowstyle='->', color='red'),
                    fontweight='bold')

plt.tight_layout()
plt.savefig('heterozygosity_analysis.png', dpi=300, bbox_inches='tight')
plt.show()

# Save detailed results
het_summary = het_data[['FID', 'IID', 'O(HOM)', 'N(NM)', 'het_rate']].copy()
het_summary['status'] = 'Normal'
het_summary.loc[het_summary['het_rate'] > threshold_high, 'status'] = 'High_Outlier'
het_summary.loc[het_summary['het_rate'] < threshold_low, 'status'] = 'Low_Outlier'

het_summary.to_csv('heterozygosity_summary.csv', index=False)
print(f"\n✓ Detailed heterozygosity summary saved to 'heterozygosity_summary.csv'")
print(f"✓ Visualization saved as 'heterozygosity_analysis.png'")
