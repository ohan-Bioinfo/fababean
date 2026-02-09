import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from matplotlib import rcParams

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

print("Creating seaborn clustermap...")

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

# Create DataFrame for seaborn
kinship_df = pd.DataFrame(matrix, index=samples, columns=samples)

# Create clustered heatmap using seaborn
plt.figure(figsize=(14, 12))
g = sns.clustermap(kinship_df, 
                   cmap='RdYlBu_r',
                   center=0,
                   vmin=-0.5, vmax=0.5,
                   figsize=(14, 12),
                   dendrogram_ratio=0.15,
                   cbar_pos=(0.02, 0.8, 0.03, 0.18),
                   tree_kws={'linewidth': 1.5})

# Customize the plot
g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xticklabels(), 
                             rotation=45, ha='right', fontsize=10)
g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_yticklabels(), 
                             fontsize=10)
g.ax_heatmap.set_xlabel('Accession', fontsize=12, labelpad=10)
g.ax_heatmap.set_ylabel('Accession', fontsize=12, labelpad=10)

# Add title
plt.suptitle('Clustered Kinship Heatmap - Faba Bean 21 Accessions', 
             fontsize=16, fontweight='bold', y=0.95)

# Adjust colorbar label
g.cax.set_ylabel('Kinship Coefficient', rotation=270, 
                 labelpad=20, fontsize=11, fontweight='bold')

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Seaborn_Clustered.png', 
            dpi=300, bbox_inches='tight')
plt.show()

print("✓ Seaborn clustered heatmap saved as 'Diversity/Kinship/Faba_Kinship_Seaborn_Clustered.png'")

# Get clustering order from seaborn
row_order = g.dendrogram_row.reordered_ind
col_order = g.dendrogram_col.reordered_ind

print(f"\n🌳 CLUSTERING ORDER (from top to bottom):")
for i, idx in enumerate(row_order):
    print(f"  {i+1:2d}. {samples[idx]}")
