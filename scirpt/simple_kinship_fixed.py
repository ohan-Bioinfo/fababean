import numpy as np
import matplotlib.pyplot as plt

print("=== SIMPLE KINSHIP FIX ===")

# Read samples
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    all_samples = [line.split()[1] for line in f]

# Read kinship data
data = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        data.append([float(x) for x in line.split()])

# Use the number of kinship lines as the sample count
n = len(data)
samples = all_samples[:n]  # Use first n samples

print(f"Using {n} samples from kinship data")
print(f"Samples: {samples}")

# Build matrix - assuming it's lower triangular including diagonal?
matrix = np.zeros((n, n))
for i in range(n):
    for j in range(len(data[i])):
        matrix[i, j] = data[i][j]
        if i != j:  # Avoid overwriting diagonal
            matrix[j, i] = data[i][j]

print(f"Matrix shape: {matrix.shape}")
print(f"Kinship range: {matrix.min():.3f} to {matrix.max():.3f}")

# Create heatmap
plt.figure(figsize=(10, 8))
plt.imshow(matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)
plt.xticks(range(n), samples, rotation=45, ha='right')
plt.yticks(range(n), samples)
plt.xlabel('Accession')
plt.ylabel('Accession')
plt.title(f'Kinship Matrix - {n} Accessions')
plt.colorbar(label='Kinship Coefficient')
plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Simple_Fixed.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Simple fixed kinship heatmap saved")
