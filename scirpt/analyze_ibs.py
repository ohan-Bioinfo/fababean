import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Set style
plt.style.use('default')
plt.rcParams['font.family'] = 'DejaVu Sans'

# Read IBS matrix (lower triangular format)
print("Reading IBS matrix...")
with open('Diversity/IBS/Faba_IBS.mdis.id', 'r') as f:
    samples = [line.strip().split()[1] for line in f.readlines()]

n_samples = len(samples)
print(f"Found {n_samples} samples")

# Read the IBS distance matrix
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs.id', dtype=float)

# Convert lower triangular to full matrix
full_matrix = np.zeros((n_samples, n_samples))
full_matrix[np.triu_indices(n_samples, k=1)] = ibs_matrix
full_matrix = full_matrix + full_matrix.T  # Make symmetric

# Convert distance to similarity (IBS proportion)
# PLINK's distance is 1 - IBS proportion, so IBS proportion = 1 - distance
ibs_similarity = 1 - full_matrix

# Set diagonal to 1 (self-similarity)
np.fill_diagonal(ibs_similarity, 1.0)

print(f"IBS similarity range: {ibs_similarity.min():.4f} - {ibs_similarity.max():.4f}")
print(f"Mean IBS similarity: {ibs_similarity.mean():.4f}")

# Create multiple visualizations
fig = plt.figure(figsize=(18, 12))

# Plot 1: Simple IBS Heatmap
ax1 = plt.subplot2grid((2, 3), (0, 0), colspan=2)
im1 = ax1.imshow(ibs_similarity, cmap='RdYlBu_r', aspect='auto', 
                vmin=0.5, vmax=1.0)  # Typical IBS range for unrelated individuals

ax1.set_xticks(range(n_samples))
ax1.set_yticks(range(n_samples))
ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=8)
ax1.set_yticklabels(samples, fontsize=8)
ax1.set_title('IBS Similarity Heatmap\n(Identity By State)', 
              fontsize=14, fontweight='bold')
cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
cbar1.set_label('IBS Proportion', rotation=270, labelpad=15)

# Plot 2: Clustered IBS Heatmap
ax2 = plt.subplot2grid((2, 3), (0, 2))

# Convert to distance for clustering
ibs_distance = 1 - ibs_similarity
linkage_matrix = linkage(squareform(ibs_distance), method='average')

# Create dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='right', 
                    ax=ax2, leaf_font_size=8)
ax2.set_title('Sample Clustering\n(Dendrogram)', fontsize=12, fontweight='bold')

# Reorder matrix based on clustering
reorder_indices = dendro['leaves']
reordered_similarity = ibs_similarity[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

# Plot 3: Clustered Heatmap
ax3 = plt.subplot2grid((2, 3), (1, 0), colspan=3)
im3 = ax3.imshow(reordered_similarity, cmap='RdYlBu_r', aspect='auto',
                vmin=0.5, vmax=1.0)

ax3.set_xticks(range(n_samples))
ax3.set_yticks(range(n_samples))
ax3.set_xticklabels(reordered_samples, rotation=45, ha='right', fontsize=8)
ax3.set_yticklabels(reordered_samples, fontsize=8)
ax3.set_title('Clustered IBS Similarity Heatmap', 
              fontsize=14, fontweight='bold')
cbar3 = plt.colorbar(im3, ax=ax3, shrink=0.8)
cbar3.set_label('IBS Proportion', rotation=270, labelpad=15)

# Add grid for better readability
for ax in [ax1, ax3]:
    ax.set_xticks(np.arange(-0.5, n_samples, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_samples, 1), minor=True)
    ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
    ax.tick_params(which="minor", size=0)

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Comprehensive.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Comprehensive IBS analysis saved as 'Diversity/IBS/Faba_IBS_Comprehensive.png'")

# Create individual detailed heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(ibs_similarity, cmap='RdYlBu_r', aspect='auto', 
               vmin=0.5, vmax=1.0)

plt.xticks(range(n_samples), samples, rotation=45, ha='right', fontsize=10)
plt.yticks(range(n_samples), samples, fontsize=10)
plt.title('IBS Similarity Matrix - Faba Bean 21 Accessions\n(Identity By State)', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('IBS Proportion (Allele Sharing)', rotation=270, labelpad=20)

# Add values for high similarity pairs (>0.95)
for i in range(n_samples):
    for j in range(n_samples):
        if i != j and ibs_similarity[i, j] > 0.95:
            plt.text(j, i, f'{ibs_similarity[i, j]:.3f}', 
                    ha='center', va='center', fontsize=7, 
                    fontweight='bold', color='black')

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Detailed_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Detailed IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Detailed_Heatmap.png'")

# Print IBS statistics
print(f"\n📊 IBS SIMILARITY STATISTICS:")
print(f"Range: {ibs_similarity.min():.4f} - {ibs_similarity.max():.4f}")
print(f"Mean: {ibs_similarity.mean():.4f} ± {ibs_similarity.std():.4f}")

# Find most similar and least similar pairs
similar_pairs = []
for i in range(n_samples):
    for j in range(i+1, n_samples):
        similar_pairs.append((samples[i], samples[j], ibs_similarity[i, j]))

similar_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n�� MOST SIMILAR PAIRS (Top 5):")
for i in range(min(5, len(similar_pairs))):
    print(f"  {similar_pairs[i][0]} - {similar_pairs[i][1]}: {similar_pairs[i][2]:.4f}")

print(f"\n🔗 LEAST SIMILAR PAIRS (Bottom 5):")
for i in range(min(5, len(similar_pairs))):
    idx = len(similar_pairs) - 1 - i
    print(f"  {similar_pairs[idx][0]} - {similar_pairs[idx][1]}: {similar_pairs[idx][2]:.4f}")

# Save IBS matrix for future use
ibs_df = pd.DataFrame(ibs_similarity, index=samples, columns=samples)
ibs_df.to_csv('Diversity/IBS/Faba_IBS_Matrix.csv')
print(f"\n✓ IBS matrix saved as 'Diversity/IBS/Faba_IBS_Matrix.csv'")
