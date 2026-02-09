import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
import networkx as nx
from scipy.cluster import hierarchy

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read IBD data
ibd_file = "Diversity/IBD/Faba_IBD.genome"
df = pd.read_csv(ibd_file, sep='\s+')

print("=== IBD ANALYSIS SUMMARY ===")
print(f"Total pairwise comparisons: {len(df)}")
print(f"Mean PI_HAT (relatedness): {df['PI_HAT'].mean():.4f}")
print(f"Standard deviation PI_HAT: {df['PI_HAT'].std():.4f}")

# Define relatedness categories
def get_relatedness_category(pi_hat):
    if pi_hat < 0.05:
        return 'Unrelated'
    elif pi_hat < 0.125:
        return '3rd Degree'
    elif pi_hat < 0.25:
        return '2nd Degree'
    elif pi_hat < 0.375:
        return '1st Degree'
    else:
        return 'Duplicate/MZ Twin'

df['relatedness'] = df['PI_HAT'].apply(get_relatedness_category)

# Count relatedness categories
relatedness_counts = df['relatedness'].value_counts()
print("\n=== RELATEDNESS DISTRIBUTION ===")
for category, count in relatedness_counts.items():
    percentage = (count / len(df)) * 100
    print(f"{category}: {count} pairs ({percentage:.1f}%)")

# Create visualizations
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))

# Plot 1: Distribution of PI_HAT values
ax1.hist(df['PI_HAT'], bins=50, color='skyblue', edgecolor='black', alpha=0.7)
ax1.axvline(x=0.125, color='red', linestyle='--', label='2nd Degree (0.125)')
ax1.axvline(x=0.25, color='orange', linestyle='--', label='1st Degree (0.25)')
ax1.axvline(x=0.375, color='purple', linestyle='--', label='Duplicate (0.375)')
ax1.set_title('Distribution of PI_HAT Values\n(Identity By Descent)', 
              fontsize=14, fontweight='bold')
ax1.set_xlabel('PI_HAT (Proportion IBD)', fontsize=12)
ax1.set_ylabel('Number of Pairs', fontsize=12)
ax1.legend()
ax1.grid(True, alpha=0.3)

# Plot 2: Relatedness categories
colors = ['#2ecc71', '#3498db', '#f39c12', '#e74c3c', '#9b59b6']
bars = ax2.bar(relatedness_counts.index, relatedness_counts.values, 
               color=colors, alpha=0.8, edgecolor='black')
ax2.set_title('Relatedness Categories Distribution', 
              fontsize=14, fontweight='bold')
ax2.set_ylabel('Number of Pairs', fontsize=12)
ax2.tick_params(axis='x', rotation=45)

# Add values on bars
for bar, count in zip(bars, relatedness_counts.values):
    height = bar.get_height()
    ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
             f'{count}', ha='center', va='bottom', fontweight='bold')

# Plot 3: Z0 vs Z1 plot (IBD segments)
scatter = ax3.scatter(df['Z0'], df['Z1'], c=df['PI_HAT'], 
                     cmap='viridis', alpha=0.6, s=50)
ax3.set_title('IBD Sharing: Z0 vs Z1\n(0 vs 1 allele IBD)', 
              fontsize=14, fontweight='bold')
ax3.set_xlabel('Z0 (Prob. sharing 0 alleles IBD)', fontsize=12)
ax3.set_ylabel('Z1 (Prob. sharing 1 allele IBD)', fontsize=12)
plt.colorbar(scatter, ax=ax3, label='PI_HAT')

# Add theoretical lines
x = np.linspace(0, 1, 100)
ax3.plot(x, 1-x, 'r--', alpha=0.5, label='Z0 + Z1 = 1')
ax3.legend()

# Plot 4: Create a relatedness heatmap
# Create sample list
samples = sorted(set(df['IID1']).union(set(df['IID2'])))
n_samples = len(samples)

# Create relatedness matrix
relatedness_matrix = np.zeros((n_samples, n_samples))
sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}

# Fill the matrix
for _, row in df.iterrows():
    i = sample_to_idx[row['IID1']]
    j = sample_to_idx[row['IID2']]
    relatedness_matrix[i, j] = row['PI_HAT']
    relatedness_matrix[j, i] = row['PI_HAT']

# Create heatmap
im = ax4.imshow(relatedness_matrix, cmap='YlOrRd', aspect='auto')
ax4.set_title('Pairwise Relatedness Heatmap\n(PI_HAT values)', 
              fontsize=14, fontweight='bold')
ax4.set_xlabel('Samples', fontsize=12)
ax4.set_ylabel('Samples', fontsize=12)

# Set ticks for better readability (show every 2nd sample if many)
if n_samples > 10:
    tick_interval = max(1, n_samples // 10)
    ticks = range(0, n_samples, tick_interval)
    ax4.set_xticks(ticks)
    ax4.set_yticks(ticks)
    ax4.set_xticklabels([samples[i] for i in ticks], rotation=45)
    ax4.set_yticklabels([samples[i] for i in ticks])
else:
    ax4.set_xticks(range(n_samples))
    ax4.set_yticks(range(n_samples))
    ax4.set_xticklabels(samples, rotation=45)
    ax4.set_yticklabels(samples)

plt.colorbar(im, ax=ax4, label='PI_HAT')

plt.tight_layout()
plt.savefig('Diversity/IBD/Faba_IBD_Analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ IBD analysis visualization saved as 'Diversity/IBD/Faba_IBD_Analysis.png'")

# Save summary statistics
summary_stats = {
    'total_pairs': len(df),
    'mean_pi_hat': df['PI_HAT'].mean(),
    'std_pi_hat': df['PI_HAT'].std(),
    'min_pi_hat': df['PI_HAT'].min(),
    'max_pi_hat': df['PI_HAT'].max(),
    'unrelated_pairs': len(df[df['PI_HAT'] < 0.05]),
    'related_pairs': len(df[df['PI_HAT'] >= 0.05])
}

with open('Diversity/IBD/IBD_summary_stats.txt', 'w') as f:
    f.write("IBD Analysis Summary Statistics\n")
    f.write("================================\n")
    for key, value in summary_stats.items():
        f.write(f"{key}: {value}\n")
    
    f.write("\nRelatedness Categories:\n")
    for category, count in relatedness_counts.items():
        percentage = (count / len(df)) * 100
        f.write(f"{category}: {count} ({percentage:.1f}%)\n")

print("✓ IBD summary statistics saved as 'Diversity/IBD/IBD_summary_stats.txt'")

# Identify closely related pairs
close_related = df[df['PI_HAT'] >= 0.125]
if len(close_related) > 0:
    print(f"\n⚠️  CLOSELY RELATED PAIRS (PI_HAT ≥ 0.125):")
    for _, pair in close_related.iterrows():
        print(f"  {pair['IID1']} - {pair['IID2']}: PI_HAT = {pair['PI_HAT']:.3f}")
else:
    print(f"\n✅ No closely related pairs detected (all PI_HAT < 0.125)")

# Check for potential duplicates
duplicates = df[df['PI_HAT'] >= 0.9]
if len(duplicates) > 0:
    print(f"\n🚨 POTENTIAL DUPLICATES (PI_HAT ≥ 0.9):")
    for _, pair in duplicates.iterrows():
        print(f"  {pair['IID1']} - {pair['IID2']}: PI_HAT = {pair['PI_HAT']:.3f}")
else:
    print(f"\n✅ No potential duplicates detected")
