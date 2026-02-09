import numpy as np
import matplotlib.pyplot as plt

# Read data
samples = [line.strip().split()[1] for line in open('Diversity/IBS/Faba_IBS.mibs.id')]
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

# Create heatmap with Accession labels
plt.figure(figsize=(12, 10))
im = plt.imshow(ibs_matrix, cmap='RdYlBu_r', vmin=0.75, vmax=1.0)

# Set ticks with "Accession" labels
plt.xticks(range(len(samples)), [f"Accession {s}" for s in samples], rotation=45, ha='right', fontsize=10)
plt.yticks(range(len(samples)), [f"Accession {s}" for s in samples], fontsize=10)

plt.title('IBS Similarity Heatmap - Faba Bean 21 Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('IBS Proportion', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Accession_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ IBS heatmap with Accession labels saved as 'Diversity/IBS/Faba_IBS_Accession_Heatmap.png'")

# Print quick stats
print(f"📊 IBS Statistics:")
print(f"Range: {ibs_matrix.min():.3f} - {ibs_matrix.max():.3f}")
print(f"Mean: {ibs_matrix.mean():.3f}")
print(f"Standard deviation: {ibs_matrix.std():.3f}")
