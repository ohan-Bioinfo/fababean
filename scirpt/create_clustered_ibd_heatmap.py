import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Read IBD data
df = pd.read_csv("Diversity/IBD/Faba_IBD.genome", sep='\s+')

# Create matrix
samples = sorted(set(df['IID1']).union(set(df['IID2'])))
n_samples = len(samples)
sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}

matrix = np.zeros((n_samples, n_samples))
for _, row in df.iterrows():
    i = sample_to_idx[row['IID1']]
    j = sample_to_idx[row['IID2']]
    matrix[i, j] = matrix[j, i] = row['PI_HAT']
np.fill_diagonal(matrix, 1.0)

# Convert to distance matrix (1 - PI_HAT) for clustering
distance_matrix = 1 - matrix

# Perform hierarchical clustering
linkage_matrix = linkage(squareform(distance_matrix), method='average')

# Create clustered heatmap
fig = plt.figure(figsize=(14, 10))

# Create gridspec
gs = plt.GridSpec(1, 2, width_ratios=[1, 4], wspace=0.01)

# Plot 1: Dendrogram
ax1 = plt.subplot(gs[0])
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax1, leaf_font_size=10)
ax1.set_title('Clustering', fontsize=12, fontweight='bold')
plt.setp(ax1.get_xticklabels(), visible=False)
ax1.spines['top'].set_visible(False)
ax1.spines['right'].set_visible(False)
ax1.spines['bottom'].set_visible(False)

# Get reordered indices from dendrogram
reorder_indices = dendro['leaves']
reordered_matrix = matrix[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

# Plot 2: Heatmap
ax2 = plt.subplot(gs[1])
im = ax2.imshow(reordered_matrix, cmap='YlOrRd', aspect='auto', vmin=0, vmax=1)

# Set ticks
ax2.set_xticks(range(n_samples))
ax2.set_yticks(range(n_samples))
ax2.set_xticklabels(reordered_samples, rotation=45, ha='right')
ax2.set_yticklabels(reordered_samples)

ax2.set_title('IBD Relatedness (Clustered)', fontsize=14, fontweight='bold')

# Add colorbar
cbar = plt.colorbar(im, ax=ax2, shrink=0.8)
cbar.set_label('PI_HAT', rotation=270, labelpad=15)

plt.suptitle('Clustered IBD Heatmap - Faba Bean 21 Accessions', 
             fontsize=16, fontweight='bold', y=0.95)

plt.tight_layout()
plt.savefig('Diversity/IBD/Faba_IBD_Clustered_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Clustered IBD heatmap saved as 'Diversity/IBD/Faba_IBD_Clustered_Heatmap.png'")

# Print cluster information
print(f"\n🌳 CLUSTERING ANALYSIS:")
print("Samples grouped by genetic similarity (IBD)")
print("Clustering method: Average linkage")
print("Distance metric: 1 - PI_HAT")
