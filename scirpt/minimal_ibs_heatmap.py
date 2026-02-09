import numpy as np
import matplotlib.pyplot as plt

# Read data
samples = [line.strip().split()[1] for line in open('Diversity/IBS/Faba_IBS.mibs.id')]
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

# Create heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(ibs_matrix, cmap='RdYlBu_r', vmin=0.75, vmax=1.0)

plt.xticks(range(len(samples)), samples, rotation=45, ha='right')
plt.yticks(range(len(samples)), samples)
plt.title('IBS Similarity Heatmap - Faba Bean 21 Accessions', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('IBS Proportion', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Minimal_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Minimal IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Minimal_Heatmap.png'")

# Print quick stats
print(f"IBS Range: {ibs_matrix.min():.3f} - {ibs_matrix.max():.3f}")
print(f"Mean IBS: {ibs_matrix.mean():.3f}")
