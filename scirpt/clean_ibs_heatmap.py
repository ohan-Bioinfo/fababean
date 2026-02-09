import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10

# Read data
samples = [line.strip().split()[1] for line in open('Diversity/IBS/Faba_IBS.mibs.id')]
ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

# Create figure
fig, ax = plt.subplots(figsize=(13, 10))

# Create heatmap with better color map
im = ax.imshow(ibs_matrix, cmap='viridis', vmin=0.75, vmax=1.0, aspect='auto')

# Set ticks with sample IDs
ax.set_xticks(range(len(samples)))
ax.set_yticks(range(len(samples)))
ax.set_xticklabels(samples, rotation=45, ha='right', fontsize=9)
ax.set_yticklabels(samples, fontsize=9)

# Set axis labels to "Accession"
ax.set_xlabel('Accession', fontsize=12, fontweight='bold', labelpad=12)
ax.set_ylabel('Accession', fontsize=12, fontweight='bold', labelpad=12)

# Set title
ax.set_title('Identity By State (IBS) Similarity\nFaba Bean Accessions', 
             fontsize=16, fontweight='bold', pad=25)

# Add colorbar
cbar = plt.colorbar(im, ax=ax, shrink=0.8, pad=0.02)
cbar.set_label('IBS Proportion', rotation=270, labelpad=20, fontsize=11, fontweight='bold')

# Add value annotations for high similarity (>0.95)
for i in range(len(samples)):
    for j in range(len(samples)):
        if i != j and ibs_matrix[i, j] > 0.95:
            ax.text(j, i, f'{ibs_matrix[i, j]:.2f}', 
                   ha='center', va='center', fontsize=6, 
                   fontweight='bold', color='white')

plt.tight_layout()
plt.savefig('Diversity/IBS/Faba_IBS_Clean_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Clean IBS heatmap saved as 'Diversity/IBS/Faba_IBS_Clean_Heatmap.png'")
