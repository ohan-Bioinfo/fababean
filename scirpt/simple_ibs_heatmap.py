import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram

# Set style
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'

# Read sample names
with open('Diversity/IBS/Faba_IBS.mibs.id', 'r') as f:
    samples = [line.strip().split()[1] for line in f.readlines()]

n_samples = len(samples)
print(f"Found {n_samples} samples")

# Read the IBS matrix - it's already a full square matrix
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

print(f"IBS matrix shape: {ibs_matrix.shape}")
print(f"IBS range: {ibs_matrix.min():.4f} - {ibs_matrix.max():.4f}")
print(f"Mean IBS: {ibs_matrix.mean():.4f}")

# Create the heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Plot 1: Simple heatmap
im1 = ax1.imshow(ibs_matrix, cmap='RdYlBu_r', aspect='auto', 
                vmin=0.75, vmax=1.0)
ax1.set_xticks(range(n_samples))
ax1.set_yticks(range(n_samples))
ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
ax1.set_yticklabels(samples, fontsize=10)
ax1.set_title('IBS Similarity Heatmap\nFaba Bean 21 Accessions', 
              fontsize=14, fontweight='bold')
cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
cbar1.set_label('IBS Proportion', rotation=270, labelpad=15)

# Add grid
ax1.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
ax1.set_yticks(np.arange(-0.5, n_samples, 1), minor=True)
ax1.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
ax1.tick_params(which="minor", size=0)

# Plot 2: Clustered heatmap
# Convert to distance for clustering
ibs_distance = 1 - ibs_matrix
linkage_matrix = linkage(ibs_distance, method='average')

# Create dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax2, leaf_font_size=10)
ax2.set_title('Sample Clustering', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Simple_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Simple IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Simple_Heatmap.png'")

# Create detailed clustered heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), 
                               gridspec_kw={'width_ratios': [1, 4]})

# Dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax1, leaf_font_size=10)
ax1.set_title('Clustering', fontsize=12, fontweight='bold')

# Reorder matrix based on clustering
reorder_indices = dendro['leaves']
reordered_matrix = ibs_matrix[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

# Clustered heatmap
im2 = ax2.imshow(reordered_matrix, cmap='RdYlBu_r', aspect='auto',
                vmin=0.75, vmax=1.0)
ax2.set_xticks(range(n_samples))
ax2.set_yticks(range(n_samples))
ax2.set_xticklabels(reordered_samples, rotation=45, ha='right', fontsize=10)
ax2.set_yticklabels(reordered_samples, fontsize=10)
ax2.set_title('Clustered IBS Similarity', fontsize=14, fontweight='bold')
cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
cbar2.set_label('IBS Proportion', rotation=270, labelpad=15)

# Add grid
ax2.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
ax2.set_yticks(np.arange(-0.5, n_samples, 1), minor=True)
ax2.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
ax2.tick_params(which="minor", size=0)

plt.suptitle('IBS Similarity - Faba Bean 21 Accessions', 
             fontsize=16, fontweight='bold', y=0.95)
plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Clustered_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Clustered IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Clustered_Heatmap.png'")

# Print statistics
print(f"\n📊 IBS STATISTICS:")
print(f"Range: {ibs_matrix.min():.4f} - {ibs_matrix.max():.4f}")
print(f"Mean: {ibs_matrix.mean():.4f} ± {ibs_matrix.std():.4f}")

# Find most similar pairs (excluding diagonal)
similar_pairs = []
for i in range(n_samples):
    for j in range(i+1, n_samples):
        similar_pairs.append((samples[i], samples[j], ibs_matrix[i, j]))

similar_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n🔗 MOST SIMILAR PAIRS (Top 5):")
for i in range(min(5, len(similar_pairs))):
    print(f"  {similar_pairs[i][0]} - {similar_pairs[i][1]}: {similar_pairs[i][2]:.4f}")

# Check for potential duplicates
print(f"\n🔍 POTENTIAL DUPLICATES (IBS > 0.99):")
duplicates = [(s1, s2, sim) for s1, s2, sim in similar_pairs if sim > 0.99]
if duplicates:
    for s1, s2, sim in duplicates:
        print(f"  {s1} - {s2}: {sim:.4f}")
else:
    print("  None found")

# Save the IBS matrix
ibs_df = pd.DataFrame(ibs_matrix, index=samples, columns=samples)
ibs_df.to_csv('Diversity/IBS/Faba_IBS_Matrix.csv')
print(f"\n✓ IBS matrix saved as 'Diversity/IBS/Faba_IBS_Matrix.csv'")
