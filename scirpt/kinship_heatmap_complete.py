import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

# Read kinship matrix (lower triangular format)
print("Reading kinship matrix...")

# Read sample IDs
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.strip().split()[1] for line in f.readlines()]

n_samples = len(samples)
print(f"Found {n_samples} accessions")

# Read kinship matrix (lower triangular)
kinship_data = np.loadtxt('Diversity/Kinship/Faba_Kinship.king', dtype=float)

# Convert lower triangular to full matrix
full_matrix = np.zeros((n_samples, n_samples))
idx = 0
for i in range(n_samples):
    for j in range(i):
        full_matrix[i, j] = kinship_data[idx]
        full_matrix[j, i] = kinship_data[idx]
        idx += 1

# Set diagonal to 1 (self-kinship)
np.fill_diagonal(full_matrix, 1.0)

print(f"Kinship range: {full_matrix.min():.4f} - {full_matrix.max():.4f}")
print(f"Mean kinship: {full_matrix.mean():.4f}")

# Create the heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Plot 1: Simple kinship heatmap
im1 = ax1.imshow(full_matrix, cmap='RdYlBu_r', aspect='auto', 
                vmin=-0.5, vmax=0.5)  # Typical kinship range

ax1.set_xticks(range(n_samples))
ax1.set_yticks(range(n_samples))
ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
ax1.set_yticklabels(samples, fontsize=10)
ax1.set_xlabel('Accession', fontsize=12, labelpad=10)
ax1.set_ylabel('Accession', fontsize=12, labelpad=10)
ax1.set_title('Kinship Matrix Heatmap\nFaba Bean Germplasm', 
              fontsize=14, fontweight='bold')

cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
cbar1.set_label('Kinship Coefficient', rotation=270, labelpad=20)

# Plot 2: Clustered kinship heatmap
# Convert to distance for clustering (1 - absolute kinship)
kinship_distance = 1 - np.abs(full_matrix)
linkage_matrix = linkage(squareform(kinship_distance), method='average')

# Create dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax2, leaf_font_size=10)
ax2.set_title('Accession Clustering\n(Dendrogram)', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Kinship analysis saved as 'Diversity/Kinship/Faba_Kinship_Analysis.png'")

# Create detailed clustered heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8), 
                               gridspec_kw={'width_ratios': [1, 4]})

# Dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax1, leaf_font_size=10)
ax1.set_title('Clustering', fontsize=12, fontweight='bold')

# Reorder matrix based on clustering
reorder_indices = dendro['leaves']
reordered_matrix = full_matrix[reorder_indices, :][:, reorder_indices]
reordered_samples = [samples[i] for i in reorder_indices]

# Clustered heatmap
im2 = ax2.imshow(reordered_matrix, cmap='RdYlBu_r', aspect='auto',
                vmin=-0.5, vmax=0.5)
ax2.set_xticks(range(n_samples))
ax2.set_yticks(range(n_samples))
ax2.set_xticklabels(reordered_samples, rotation=45, ha='right', fontsize=10)
ax2.set_yticklabels(reordered_samples, fontsize=10)
ax2.set_xlabel('Accession', fontsize=12, labelpad=10)
ax2.set_ylabel('Accession', fontsize=12, labelpad=10)
ax2.set_title('Clustered Kinship Matrix', fontsize=14, fontweight='bold')
cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
cbar2.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.suptitle('Kinship Analysis - Faba Bean 21 Accessions', 
             fontsize=16, fontweight='bold', y=0.95)
plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Clustered.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Clustered kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Clustered.png'")

# Print kinship statistics
print(f"\n📊 KINSHIP STATISTICS:")
print(f"Range: {full_matrix.min():.4f} - {full_matrix.max():.4f}")
print(f"Mean: {full_matrix.mean():.4f} ± {full_matrix.std():.4f}")

# Find most related and least related pairs
related_pairs = []
for i in range(n_samples):
    for j in range(i+1, n_samples):
        related_pairs.append((samples[i], samples[j], full_matrix[i, j]))

related_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n🔗 MOST RELATED PAIRS (Top 5):")
for i in range(min(5, len(related_pairs))):
    print(f"  {related_pairs[i][0]} - {related_pairs[i][1]}: {related_pairs[i][2]:.4f}")

print(f"\n🔗 LEAST RELATED PAIRS (Bottom 5):")
for i in range(min(5, len(related_pairs))):
    idx = len(related_pairs) - 1 - i
    print(f"  {related_pairs[idx][0]} - {related_pairs[idx][1]}: {related_pairs[idx][2]:.4f}")

# Interpret kinship values
print(f"\n💡 KINSHIP INTERPRETATION:")
print("≈ 0.5: Parent-Offspring / Full Siblings")
print("≈ 0.25: Half Siblings / Grandparent-Grandchild")
print("≈ 0.125: First Cousins")
print("≈ 0.0625: Second Cousins")
print("≈ 0: Unrelated")
print("< 0: Less related than population average")

# Save kinship matrix
kinship_df = pd.DataFrame(full_matrix, index=samples, columns=samples)
kinship_df.to_csv('Diversity/Kinship/Faba_Kinship_Matrix.csv')
print(f"\n✓ Kinship matrix saved as 'Diversity/Kinship/Faba_Kinship_Matrix.csv'")
