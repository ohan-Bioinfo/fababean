import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd
from scipy.cluster.hierarchy import linkage, dendrogram

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

print("Reading kinship matrix...")

# Read sample IDs
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.strip().split()[1] for line in f.readlines()]

n_samples = len(samples)
print(f"Found {n_samples} accessions")

# Read kinship matrix line by line (ragged format)
kinship_lines = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        kinship_lines.append([float(x) for x in line.strip().split()])

print(f"Kinship file has {len(kinship_lines)} lines")

# Build the full symmetric matrix
full_matrix = np.zeros((n_samples, n_samples))

# Fill the lower triangular part
for i in range(n_samples):
    for j in range(len(kinship_lines[i])):
        full_matrix[i, j] = kinship_lines[i][j]
        full_matrix[j, i] = kinship_lines[i][j]  # Make symmetric

# Note: The diagonal is already included in the KING output

print(f"Kinship range: {full_matrix.min():.4f} - {full_matrix.max():.4f}")
print(f"Mean kinship: {full_matrix.mean():.4f}")

# Create minimal heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(full_matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)

plt.xticks(range(len(samples)), samples, rotation=45, ha='right')
plt.yticks(range(len(samples)), samples)
plt.xlabel('Accession', fontsize=12, labelpad=10)
plt.ylabel('Accession', fontsize=12, labelpad=10)
plt.title('Kinship Matrix - Faba Bean Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Fixed_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Fixed kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Fixed_Heatmap.png'")

# Print statistics
print(f"\n📊 KINSHIP STATISTICS:")
print(f"Accessions: {n_samples}")
print(f"Kinship Range: {full_matrix.min():.3f} - {full_matrix.max():.3f}")
print(f"Mean Kinship: {full_matrix.mean():.3f}")
print(f"Standard Deviation: {full_matrix.std():.3f}")

# Find most related pairs
related_pairs = []
for i in range(n_samples):
    for j in range(i+1, n_samples):
        related_pairs.append((samples[i], samples[j], full_matrix[i, j]))

related_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n🔗 MOST RELATED PAIRS (Top 5):")
for i in range(min(5, len(related_pairs))):
    print(f"  {related_pairs[i][0]} - {related_pairs[i][1]}: {related_pairs[i][2]:.4f}")

# Save the matrix
kinship_df = pd.DataFrame(full_matrix, index=samples, columns=samples)
kinship_df.to_csv('Diversity/Kinship/Faba_Kinship_Matrix.csv')
print(f"\n✓ Kinship matrix saved as 'Diversity/Kinship/Faba_Kinship_Matrix.csv'")
