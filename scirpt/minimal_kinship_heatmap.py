import numpy as np
import matplotlib.pyplot as plt

# Read kinship data
samples = [line.strip().split()[1] for line in open('Diversity/Kinship/Faba_Kinship.king.id')]
kinship_data = np.loadtxt('Diversity/Kinship/Faba_Kinship.king')

# Reconstruct full matrix
n_samples = len(samples)
full_matrix = np.zeros((n_samples, n_samples))
idx = 0
for i in range(n_samples):
    for j in range(i):
        full_matrix[i, j] = kinship_data[idx]
        full_matrix[j, i] = kinship_data[idx]
        idx += 1
np.fill_diagonal(full_matrix, 1.0)

# Create minimal heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(full_matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)

plt.xticks(range(len(samples)), samples, rotation=45, ha='right')
plt.yticks(range(len(samples)), samples)
plt.xlabel('Accession', fontsize=12, labelpad=10)
plt.ylabel('Accession', fontsize=12, labelpad=10)
plt.title('Kinship Matrix - Faba Bean 21 Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Minimal_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Minimal kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Minimal_Heatmap.png'")

# Print summary
print(f"Kinship Range: {full_matrix.min():.3f} - {full_matrix.max():.3f}")
print(f"Mean Kinship: {full_matrix.mean():.3f}")
