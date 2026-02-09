import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10

# Read IBD data
ibd_file = "Diversity/IBD/Faba_IBD.genome"
df = pd.read_csv(ibd_file, sep='\s+')

# Get unique samples and sort them
samples = sorted(set(df['IID1']).union(set(df['IID2'])))
n_samples = len(samples)
print(f"Creating IBD heatmap for {n_samples} samples")

# Create relatedness matrix
relatedness_matrix = np.zeros((n_samples, n_samples))
sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}

# Fill the matrix with PI_HAT values
for _, row in df.iterrows():
    i = sample_to_idx[row['IID1']]
    j = sample_to_idx[row['IID2']]
    relatedness_matrix[i, j] = row['PI_HAT']
    relatedness_matrix[j, i] = row['PI_HAT']

# Set diagonal to 1 (self-relatedness)
np.fill_diagonal(relatedness_matrix, 1.0)

# Create the heatmap
fig, ax = plt.subplots(figsize=(10, 8))

# Create heatmap with custom color scale
im = ax.imshow(relatedness_matrix, cmap='YlOrRd', aspect='auto', 
               vmin=0, vmax=1, interpolation='nearest')

# Set ticks
ax.set_xticks(range(n_samples))
ax.set_yticks(range(n_samples))
ax.set_xticklabels(samples, rotation=45, ha='right')
ax.set_yticklabels(samples)

# Add title and labels
ax.set_title('Identity By Descent (IBD) Heatmap\nFaba Bean - 21 Accessions', 
             fontsize=14, fontweight='bold', pad=20)
ax.set_xlabel('Accession', fontsize=12, labelpad=10)
ax.set_ylabel('Accession', fontsize=12, labelpad=10)

# Add colorbar with label
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('PI_HAT (Proportion IBD)', fontsize=11, rotation=270, labelpad=15)

# Add grid for better readability
ax.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
ax.set_yticks(np.arange(-0.5, n_samples, 1), minor=True)
ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
ax.tick_params(which="minor", size=0)

# Add PI_HAT values on the heatmap (only for values > 0.1 for clarity)
for i in range(n_samples):
    for j in range(n_samples):
        value = relatedness_matrix[i, j]
        if value > 0.1 and i != j:  # Don't show self-relatedness (1.0)
            color = 'white' if value > 0.5 else 'black'
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', 
                    color=color, fontsize=7, fontweight='bold')

plt.tight_layout()
plt.savefig('Diversity/IBD/Faba_IBD_Heatmap_21samples.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ IBD heatmap saved as 'Diversity/IBD/Faba_IBD_Heatmap_21samples.png'")

# Print summary of relatedness
print(f"\n📊 IBD SUMMARY FOR 21 ACCESSIONS:")
print(f"Mean PI_HAT: {df['PI_HAT'].mean():.4f}")
print(f"PI_HAT range: {df['PI_HAT'].min():.4f} - {df['PI_HAT'].max():.4f}")

# Count related pairs
related_pairs = len(df[df['PI_HAT'] >= 0.125])
total_pairs = len(df)
print(f"Pairs with PI_HAT ≥ 0.125 (2nd degree+): {related_pairs}/{total_pairs} ({(related_pairs/total_pairs)*100:.1f}%)")

# Identify top related pairs
if related_pairs > 0:
    print(f"\n🔗 TOP RELATED PAIRS:")
    top_related = df[df['PI_HAT'] >= 0.125].nlargest(5, 'PI_HAT')
    for _, row in top_related.iterrows():
        print(f"  {row['IID1']} - {row['IID2']}: PI_HAT = {row['PI_HAT']:.3f}")
