import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set professional style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10

# Read data
samples = [line.strip().split()[1] for line in open('Diversity/IBS/Faba_IBS.mibs.id')]
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

# Create figure with professional layout
fig, ax = plt.subplots(figsize=(14, 11))

# Create heatmap
im = ax.imshow(ibs_matrix, cmap='RdYlBu_r', vmin=0.75, vmax=1.0, aspect='auto')

# Set labels with "Accession" prefix
accession_labels = [f"Accession {s}" for s in samples]
ax.set_xticks(range(len(samples)))
ax.set_yticks(range(len(samples)))
ax.set_xticklabels(accession_labels, rotation=45, ha='right', fontsize=9)
ax.set_yticklabels(accession_labels, fontsize=9)

# Set titles and labels
ax.set_title('Identity By State (IBS) Similarity Matrix\nFaba Bean Germplasm Collection', 
             fontsize=16, fontweight='bold', pad=25)
ax.set_xlabel('Accessions', fontsize=12, labelpad=15)
ax.set_ylabel('Accessions', fontsize=12, labelpad=15)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('IBS Proportion\n(Allele Sharing)', rotation=270, labelpad=20, fontsize=11)

# Add grid for better readability
ax.set_xticks(np.arange(-0.5, len(samples), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(samples), 1), minor=True)
ax.grid(which="minor", color="white", linestyle='-', linewidth=0.8)
ax.tick_params(which="minor", size=0)

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Professional_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Professional IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Professional_Heatmap.png'")

# Calculate and display key statistics
n = len(samples)
print(f"\n📊 ANALYSIS SUMMARY:")
print(f"Accessions: {n}")
print(f"IBS Range: {ibs_matrix.min():.4f} - {ibs_matrix.max():.4f}")
print(f"Mean Similarity: {ibs_matrix.mean():.4f}")

# Find most distinct and most similar accessions (excluding self-comparisons)
similar_pairs = []
for i in range(n):
    for j in range(i+1, n):
        similar_pairs.append((samples[i], samples[j], ibs_matrix[i, j]))

if similar_pairs:
    similar_pairs.sort(key=lambda x: x[2], reverse=True)
    print(f"\n🔗 MOST SIMILAR ACCESSIONS:")
    for i in range(min(3, len(similar_pairs))):
        print(f"  Accession {similar_pairs[i][0]} - Accession {similar_pairs[i][1]}: {similar_pairs[i][2]:.4f}")
    
    print(f"\n🔗 MOST DISTINCT ACCESSIONS:")
    for i in range(min(3, len(similar_pairs))):
        idx = len(similar_pairs) - 1 - i
        print(f"  Accession {similar_pairs[idx][0]} - Accession {similar_pairs[idx][1]}: {similar_pairs[idx][2]:.4f}")
