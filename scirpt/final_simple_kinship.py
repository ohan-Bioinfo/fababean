import numpy as np
import matplotlib.pyplot as plt

print("=== FINAL SIMPLE KINSHIP ===")

# Read samples (skip header)
with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
    samples = [line.split()[1] for line in f.readlines()[1:]]  # Skip header

print(f"Actual samples: {len(samples)}")
print(f"Sample list: {samples}")

# Read kinship data
kinship_data = []
with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
    for line in f:
        kinship_data.append([float(x) for x in line.split()])

print(f"Kinship lines: {len(kinship_data)}")

# Build matrix for 21 samples
n = 21
matrix = np.zeros((n, n))

# Fill lower triangular (without diagonal)
for i in range(len(kinship_data)):
    for j in range(len(kinship_data[i])):
        matrix[i+1, j] = kinship_data[i][j]
        matrix[j, i+1] = kinship_data[i][j]  # Symmetric

# Set diagonal to 1
np.fill_diagonal(matrix, 1.0)

print(f"Matrix shape: {matrix.shape}")
print(f"Kinship range: {matrix.min():.3f} to {matrix.max():.3f}")

# Create final heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(matrix, cmap='RdYlBu_r', vmin=-0.5, vmax=0.5)

plt.xticks(range(n), samples, rotation=45, ha='right')
plt.yticks(range(n), samples)
plt.xlabel('Accession', fontsize=12, labelpad=10)
plt.ylabel('Accession', fontsize=12, labelpad=10)
plt.title('Kinship Matrix - Faba Bean 21 Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('Kinship Coefficient', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/Kinship/Faba_Kinship_Final.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Final kinship heatmap saved as 'Diversity/Kinship/Faba_Kinship_Final.png'")

# Print summary
print(f"\n📊 FINAL KINSHIP SUMMARY:")
print(f"Accessions: {n}")
print(f"Kinship Range: {matrix.min():.3f} - {matrix.max():.3f}")
print(f"Mean: {matrix.mean():.3f}")

# Check for close relatives
close_relatives = []
for i in range(n):
    for j in range(i+1, n):
        if matrix[i,j] > 0.1:  # Threshold for relatedness
            close_relatives.append((samples[i], samples[j], matrix[i,j]))

if close_relatives:
    print(f"\n🔍 CLOSE RELATIVES (kinship > 0.1):")
    for pair in sorted(close_relatives, key=lambda x: x[2], reverse=True)[:5]:
        print(f"  {pair[0]} - {pair[1]}: {pair[2]:.3f}")
else:
    print(f"\n✅ No close relatives detected (all kinship ≤ 0.1)")
