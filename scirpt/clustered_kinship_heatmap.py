import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

print("Creating clustered kinship heatmap with dendrograms...")

# Read samples (skip header)
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.split()[1] for line in f.readlines()[1:]]

print(f"Processing {len(samples)} accessions")

# Read kinship data
kinship_data = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        kinship_data.append([float(x) for x in line.split()])

# Build kinship matrix
n = len(samples)
matrix = np.zeros((n, n))

# Fill lower triangular (without diagonal)
for i in range(len(kinship_data)):
    for j in range(len(kinship_data[i])):
        matrix[i+1, j] = kinship_data[i][j]
        matrix[j, i+1] = kinship_data[i][j]  # Symmetric

# Set diagonal to 1
np.fill_diagonal(matrix, 1.0)

print(f"Kinship matrix shape: {matrix.shape}")
print(f"Kinship range: {matrix.min():.3f} to {matrix.max():.3f}")

# Perform hierarchical clustering
# Convert to distance matrix (1 - absolute kinship)
distance_matrix = 1 - np.abs(matrix)

# Cluster rows and columns
row_linkage = linkage(squareform(distance_matrix), method='average')
col_linkage = linkage(squareform(distance_matrix.T), method='average')

# Create figure with subplots
fig = plt.figure(figsize=(16, 14))

# Define the grid layout
gs = plt.GridSpec(4, 4, figure=fig, 
                  height_ratios=[0.5, 0.1, 3, 0.5],
                  width_ratios=[0.5, 0.1, 3, 0.5],
                  hspace=0.01, wspace=0.01)

# Plot 1: Top dendrogram
ax_dendro_top = fig.add_subplot(gs[0, 2])
dendro_top = dendrogram(row_linkage, orientation='top', 
                       labels=samples, ax=ax_dendro_top,
                       color_threshold=0.7*np.max(row_linkage[:, 2]))
ax_dendro_top.set_xticks([])
ax_dendro_top.set_yticks([])
ax_dendro_top.spines['top'].set_visible(False)
ax_dendro_top.spines['right'].set_visible(False)
ax_dendro_top.spines['bottom'].set_visible(False)
ax_dendro_top.spines['left'].set_visible(False)
ax_dendro_top.set_title('Hierarchical Clustering - Faba Bean Kinship', 
                       fontsize=16, fontweight='bold', pad=20)

# Plot 2: Left dendrogram
ax_dendro_left = fig.add_subplot(gs[2, 0])
dendro_left = dendrogram(row_linkage, orientation='left', 
                        labels=samples, ax=ax_dendro_left,
                        color_threshold=0.7*np.max(row_linkage[:, 2]))
ax_dendro_left.set_xticks([])
ax_dendro_left.set_yticks([])
ax_dendro_left.spines['top'].set_visible(False)
ax_dendro_left.spines['right'].set_visible(False)
ax_dendro_left.spines['bottom'].set_visible(False)
ax_dendro_left.spines['left'].set_visible(False)

# Plot 3: Main heatmap
ax_heatmap = fig.add_subplot(gs[2, 2])

# Reorder the matrix based on clustering
reorder_indices = dendro_left['leaves']
reordered_matrix = matrix[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

# Create the heatmap
im = ax_heatmap.imshow(reordered_matrix, cmap='RdYlBu_r', 
                      aspect='auto', vmin=-0.5, vmax=0.5)

# Set ticks for the heatmap
ax_heatmap.set_xticks(range(n))
ax_heatmap.set_yticks(range(n))
ax_heatmap.set_xticklabels(reordered_samples, rotation=45, 
                          ha='right', fontsize=10)
ax_heatmap.set_yticklabels(reordered_samples, fontsize=10)
ax_heatmap.set_xlabel('Accession', fontsize=12, labelpad=10)
ax_heatmap.set_ylabel('Accession', fontsize=12, labelpad=10)

# Add grid for better readability
ax_heatmap.set_xticks(np.arange(-0.5, n, 1), minor=True)
ax_heatmap.set_yticks(np.arange(-0.5, n, 1), minor=True)
ax_heatmap.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
ax_heatmap.tick_params(which="minor", size=0)

# Plot 4: Colorbar
ax_colorbar = fig.add_subplot(gs[2, 3])
cbar = plt.colorbar(im, cax=ax_colorbar)
cbar.set_label('Kinship Coefficient', rotation=270, 
               labelpad=20, fontsize=12, fontweight='bold')

# Add some statistics text
ax_stats = fig.add_subplot(gs[3, 2])
ax_stats.axis('off')
stats_text = f'Kinship Range: {matrix.min():.3f} to {matrix.max():.3f} | Mean: {matrix.mean():.3f} | Accessions: {n}'
ax_stats.text(0.5, 0.5, stats_text, ha='center', va='center', 
              fontsize=10, transform=ax_stats.transAxes)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Clustered_Complete.png', 
            dpi=300, bbox_inches='tight')
plt.show()

print("✓ Complete clustered kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Clustered_Complete.png'")

# Print clustering information
print(f"\n📊 CLUSTERING ANALYSIS:")
print(f"Number of clusters detected: {len(set(dendro_left['color_list']))}")
print(f"Accessions in order of clustering:")
for i, idx in enumerate(reorder_indices):
    print(f"  {i+1:2d}. {samples[idx]}")
