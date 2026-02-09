import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

print("=== ROBUST KINSHIP READER ===")

# Read sample IDs
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    all_samples = [line.strip().split()[1] for line in f.readlines()]

print(f"All samples in .king.id: {len(all_samples)}")
print(f"Sample names: {all_samples}")

# Read kinship data
kinship_lines = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        if line.strip():  # Skip empty lines
            kinship_lines.append([float(x) for x in line.strip().split()])

print(f"Kinship file lines: {len(kinship_lines)}")

# Determine the actual number of samples in kinship data
# For a lower triangular matrix without diagonal, we have n*(n-1)/2 elements
# Let's find n such that n*(n-1)/2 = total_elements
total_elements = sum(len(line) for line in kinship_lines)
print(f"Total kinship values: {total_elements}")

# Solve for n: n*(n-1)/2 = total_elements
# n^2 - n - 2*total_elements = 0
n_kinship = int((1 + np.sqrt(1 + 8 * total_elements)) / 2)
print(f"Calculated samples from kinship data: {n_kinship}")

# Use only the first n_kinship samples
samples = all_samples[:n_kinship]
print(f"Using {len(samples)} samples: {samples}")

# Build the kinship matrix
kinship_matrix = np.zeros((n_kinship, n_kinship))

# Fill the matrix from the triangular data
idx = 0
for i in range(n_kinship):
    for j in range(i):  # j from 0 to i-1 (lower triangular without diagonal)
        kinship_matrix[i, j] = kinship_lines[i][j]
        kinship_matrix[j, i] = kinship_lines[i][j]  # Make symmetric

# Set diagonal to 0 (typical for kinship matrices) or we can set to 1 for self?
# For visualization, setting diagonal to 1 might be better
np.fill_diagonal(kinship_matrix, 1.0)

print(f"Kinship matrix shape: {kinship_matrix.shape}")
print(f"Kinship range: {kinship_matrix.min():.4f} - {kinship_matrix.max():.4f}")

# Create the heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(kinship_matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)

plt.xticks(range(len(samples)), samples, rotation=45, ha='right')
plt.yticks(range(len(samples)), samples)
plt.xlabel('Accession', fontsize=12, labelpad=10)
plt.ylabel('Accession', fontsize=12, labelpad=10)
plt.title(f'Kinship Matrix - Faba Bean {len(samples)} Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Robust_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Robust kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Robust_Heatmap.png'")

# Print statistics
print(f"\n📊 KINSHIP STATISTICS:")
print(f"Accessions: {len(samples)}")
print(f"Kinship Range: {kinship_matrix.min():.3f} - {kinship_matrix.max():.3f}")
print(f"Mean Kinship: {kinship_matrix.mean():.3f}")

# Find most related pairs
related_pairs = []
for i in range(len(samples)):
    for j in range(i+1, len(samples)):
        related_pairs.append((samples[i], samples[j], kinship_matrix[i, j]))

related_pairs.sort(key=lambda x: x[2], reverse=True)

print(f"\n🔗 MOST RELATED PAIRS (Top 5):")
for i in range(min(5, len(related_pairs))):
    print(f"  {related_pairs[i][0]} - {related_pairs[i][1]}: {related_pairs[i][2]:.4f}")

# Save the matrix
kinship_df = pd.DataFrame(kinship_matrix, index=samples, columns=samples)
kinship_df.to_csv('Diversity/Kinship/Faba_Kinship_Matrix.csv')
print(f"\n✓ Kinship matrix saved as 'Diversity/Kinship/Faba_Kinship_Matrix.csv'")
