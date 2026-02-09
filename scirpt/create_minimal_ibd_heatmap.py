import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read IBD data
df = pd.read_csv("Diversity/IBD/Faba_IBD.genome", sep='\s+')

# Create matrix
samples = sorted(set(df['IID1']).union(set(df['IID2'])))
n_samples = len(samples)
sample_to_idx = {sample: idx for idx, sample in enumerate(samples)}

matrix = np.zeros((n_samples, n_samples))
for _, row in df.iterrows():
    i = sample_to_idx[row['IID1']]
    j = sample_to_idx[row['IID2']]
    matrix[i, j] = matrix[j, i] = row['PI_HAT']
np.fill_diagonal(matrix, 1.0)

# Create clean heatmap
plt.figure(figsize=(12, 10))
im = plt.imshow(matrix, cmap='YlOrRd', vmin=0, vmax=1)

plt.xticks(range(n_samples), samples, rotation=45, ha='right')
plt.yticks(range(n_samples), samples)
plt.title('IBD Relatedness Heatmap - Faba Bean (21 Accessions)', 
          fontsize=16, fontweight='bold', pad=20)

cbar = plt.colorbar(im, shrink=0.8)
cbar.set_label('PI_HAT (Genetic Relatedness)', rotation=270, labelpad=20)

plt.tight_layout()
plt.savefig('Diversity/IBD/Faba_IBD_Minimal_Heatmap.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Minimal IBD heatmap saved as 'Diversity/IBD/Faba_IBD_Minimal_Heatmap.png'")
