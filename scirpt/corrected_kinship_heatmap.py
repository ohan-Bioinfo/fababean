import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import squareform

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

print("=== CORRECTED KINSHIP ANALYSIS ===")

# Read sample IDs - SKIP THE HEADER
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    lines = f.readlines()
    # Skip the first line (header) and take the rest
    samples = [line.strip().split()[1] for line in lines[1:]]

print(f"Number of actual samples: {len(samples)}")
print(f"Sample names: {samples}")

# Read kinship file line by line
kinship_lines = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        if line.strip():  # Skip empty lines
            kinship_lines.append([float(x) for x in line.strip().split()])

print(f"Kinship file lines: {len(kinship_lines)}")

# For n samples, we expect n-1 lines in the kinship file (lower triangular without diagonal)
# So kinship_lines should have len(samples) - 1 lines
expected_lines = len(samples) - 1
print(f"Expected kinship lines for {len(samples)} samples: {expected_lines}")

if len(kinship_lines) != expected_lines:
    print(f"⚠️  WARNING: Expected {expected_lines} lines but got {len(kinship_lines)}")
    print("Using available data to build matrix...")

# Build the kinship matrix
n = len(samples)
kinship_matrix = np.zeros((n, n))

# Fill the lower triangular part (without diagonal)
for i in range(len(kinship_lines)):
    for j in range(len(kinship_lines[i])):
        # i+1 because kinship_lines starts from row 1 (not row 0)
        kinship_matrix[i+1, j] = kinship_lines[i][j]
        kinship_matrix[j, i+1] = kinship_lines[i][j]  # Make symmetric

# Set diagonal to 1 (self-kinship for visualization)
np.fill_diagonal(kinship_matrix, 1.0)

print(f"Kinship matrix shape: {kinship_matrix.shape}")
print(f"Kinship range: {kinship_matrix.min():.4f} - {kinship_matrix.max():.4f}")
print(f"Mean kinship: {kinship_matrix.mean():.4f}")

# Create the heatmap
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))

# Plot 1: Simple kinship heatmap
im1 = ax1.imshow(kinship_matrix, cmap='RdYlBu_r', aspect='auto', 
                vmin=-0.5, vmax=0.5)

ax1.set_xticks(range(n))
ax1.set_yticks(range(n))
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
kinship_distance = 1 - np.abs(kinship_matrix)
linkage_matrix = linkage(squareform(kinship_distance), method='average')

# Create dendrogram
dendro = dendrogram(linkage_matrix, labels=samples, orientation='left', 
                    ax=ax2, leaf_font_size=10)
ax2.set_title('Accession Clustering\n(Dendrogram)', fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Corrected_Analysis.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Corrected kinship analysis saved as 'Diversity/Kinship/Faba_Kinship_Corrected_Analysis.png'")

# Create minimal heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(kinship_matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)

plt.xticks(range(n), samples, rotation=45, ha='right')
plt.yticks(range(n), samples)
plt.xlabel('Accession', fontsize=12, labelpad=10)
plt.ylabel('Accession', fontsize=12, labelpad=10)
plt.title('Kinship Matrix - Faba Bean 21 Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Corrected_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Corrected kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Corrected_Heatmap.png'")

# Print statistics
print(f"\n📊 KINSHIP STATISTICS:")
print(f"Accessions: {n}")
print(f"Kinship Range: {kinship_matrix.min():.3f} - {kinship_matrix.max():.3f}")
print(f"Mean Kinship: {kinship_matrix.mean():.3f}")
print(f"Standard Deviation: {kinship_matrix.std():.3f}")

# Find most related pairs
related_pairs = []
for i in range(n):
    for j in range(i+1, n):
        related_pairs.append((samples[i], samples[j], kinship_matrix[i, j]))

related_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n🔗 MOST RELATED PAIRS (Top 5):")
for i in range(min(5, len(related_pairs))):
    print(f"  {related_pairs[i][0]} - {related_pairs[i][1]}: {related_pairs[i][2]:.4f}")

print(f"\n🔗 LEAST RELATED PAIRS (Bottom 5):")
for i in range(min(5, len(related_pairs))):
    idx = len(related_pairs) - 1 - i
    print(f"  {related_pairs[idx][0]} - {related_pairs[idx][1]}: {related_pairs[idx][2]:.4f}")

# Save the matrix
kinship_df = pd.DataFrame(kinship_matrix, index=samples, columns=samples)
kinship_df.to_csv('Diversity/Kinship/Faba_Kinship_Corrected_Matrix.csv')
print(f"\n✓ Kinship matrix saved as 'Diversity/Kinship/Faba_Kinship_Corrected_Matrix.csv'")
