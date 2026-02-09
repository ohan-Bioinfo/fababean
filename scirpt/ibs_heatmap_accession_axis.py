import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

# Read data
samples = [line.strip().split()[1] for line in open('Diversity/IBS/Faba_IBS.mibs.id')]
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

# Create figure
fig, ax = plt.subplots(figsize=(12, 10))

# Create heatmap
im = ax.imshow(ibs_matrix, cmap='RdYlBu_r', vmin=0.75, vmax=1.0, aspect='auto')

# Set ticks with original sample IDs (F1, F2, etc.)
ax.set_xticks(range(len(samples)))
ax.set_yticks(range(len(samples)))
ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
ax.set_yticklabels(samples, fontsize=10)

# Set axis labels to "Accession" instead of "Samples"
ax.set_xlabel('Accession', fontsize=12, labelpad=10)
ax.set_ylabel('Accession', fontsize=12, labelpad=10)

# Set title
ax.set_title('IBS Similarity Heatmap\nFaba Bean Germplasm', 
             fontsize=16, fontweight='bold', pad=20)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8)
cbar.set_label('IBS Proportion', rotation=270, labelpad=20, fontsize=11)

# Add grid for better readability
ax.set_xticks(np.arange(-0.5, len(samples), 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(samples), 1), minor=True)
ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
ax.tick_params(which="minor", size=0)

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Heatmap_Accession.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ IBS heatmap with Accession axis labels saved as 'Diversity/IBS/Faba_IBS_Heatmap_Accession.png'")

# Print statistics
print(f"📊 IBS Statistics:")
print(f"Accessions: {len(samples)}")
print(f"IBS Range: {ibs_matrix.min():.3f} - {ibs_matrix.max():.3f}")
print(f"Mean Similarity: {ibs_matrix.mean():.3f}")
