import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Set publication quality style
plt.style.use('default')
rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 14,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'legend.fontsize': 10,
    'figure.titlesize': 16
})

print("Creating publication-ready clustered heatmap...")

# Read samples (skip header)
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.split()[1] for line in f.readlines()[1:]]

# Read kinship data
kinship_data = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        kinship_data.append([float(x) for x in line.split()])

# Build kinship matrix
n = len(samples)
matrix = np.zeros((n, n))
for i in range(len(kinship_data)):
    for j in range(len(kinship_data[i])):
        matrix[i+1, j] = kinship_data[i][j]
        matrix[j, i+1] = kinship_data[i][j]
np.fill_diagonal(matrix, 1.0)

# Perform clustering
distance_matrix = 1 - np.abs(matrix)
row_linkage = linkage(squareform(distance_matrix), method='average')

# Create figure with specific layout
fig = plt.figure(figsize=(10, 8))

# Create gridspec
gs = plt.GridSpec(2, 2, figure=fig, 
                  width_ratios=[0.2, 1],
                  height_ratios=[0.2, 1],
                  wspace=0.02, hspace=0.02)

# Top dendrogram
ax_dend_top = fig.add_subplot(gs[0, 1])
dend_top = dendrogram(row_linkage, orientation='top', 
                     color_threshold=0.7*np.max(row_linkage[:, 2]),
                     ax=ax_dend_top, no_labels=True)
ax_dend_top.axis('off')

# Left dendrogram
ax_dend_left = fig.add_subplot(gs[1, 0])
dend_left = dendrogram(row_linkage, orientation='left', 
                      color_threshold=0.7*np.max(row_linkage[:, 2]),
                      ax=ax_dend_left, no_labels=True)
ax_dend_left.axis('off')

# Main heatmap
ax_heat = fig.add_subplot(gs[1, 1])
reorder_indices = dend_left['leaves']
reordered_matrix = matrix[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

im = ax_heat.imshow(reordered_matrix, cmap='coolwarm', 
                   aspect='auto', vmin=-0.5, vmax=0.5)

# Set ticks
ax_heat.set_xticks(range(n))
ax_heat.set_yticks(range(n))
ax_heat.set_xticklabels(reordered_samples, rotation=45, 
                       ha='right', fontsize=8)
ax_heat.set_yticklabels(reordered_samples, fontsize=8)
ax_heat.set_xlabel('Accession', fontsize=11, fontweight='bold')
ax_heat.set_ylabel('Accession', fontsize=11, fontweight='bold')

# Add colorbar
cbar = plt.colorbar(im, ax=ax_heat, shrink=0.8, pad=0.02)
cbar.set_label('Kinship Coefficient', rotation=270, 
               labelpad=15, fontsize=10, fontweight='bold')

# Add title
plt.suptitle('Genetic Kinship Among Faba Bean Accessions', 
             fontsize=16, fontweight='bold', y=0.95)

# Add statistics box
stats_text = f'n = {n} accessions\\nRange: {matrix.min():.3f}–{matrix.max():.3f}\\nMean: {matrix.mean():.3f}'
ax_heat.text(0.02, 0.98, stats_text, transform=ax_heat.transAxes,
            fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Publication_Ready.png', 
            dpi=600, bbox_inches='tight', 
            facecolor='white', edgecolor='none')
plt.show()

print("✓ Publication-ready clustered heatmap saved as 'Diversity/Kinship/Faba_Kinship_Publication_Ready.png'")

print(f"\n📋 PUBLICATION DETAILS:")
print(f"Figure size: 10x8 inches")
print(f"Resolution: 600 DPI")
print(f"Color scheme: Coolwarm")
print(f"Font: Arial")
print(f"Clustering method: Average linkage")
