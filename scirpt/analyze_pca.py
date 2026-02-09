import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from sklearn.preprocessing import StandardScaler

# Set style
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'

# Read PCA data
eigenval_file = "Diversity/PCA/Faba_PCA.eigenval"
eigenvec_file = "Diversity/PCA/Faba_PCA.eigenvec"

# Read eigenvalues
eigenvals = pd.read_csv(eigenval_file, header=None, names=['eigenvalue'])
eigenvals['variance_explained'] = eigenvals['eigenvalue'] / eigenvals['eigenvalue'].sum() * 100
eigenvals['pc'] = [f'PC{i+1}' for i in range(len(eigenvals))]

# Read eigenvectors
pca_df = pd.read_csv(eigenvec_file, sep='\s+', header=None)
pca_df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, pca_df.shape[1]-1)]

print("=== PCA ANALYSIS SUMMARY ===")
print(f"Total variance explained by top 10 PCs: {eigenvals['variance_explained'].sum():.2f}%")
print(f"Variance explained by PC1: {eigenvals['variance_explained'].iloc[0]:.2f}%")
print(f"Variance explained by PC2: {eigenvals['variance_explained'].iloc[1]:.2f}%")
print(f"Variance explained by PC3: {eigenvals['variance_explained'].iloc[2]:.2f}%")

# Create visualizations
fig = plt.figure(figsize=(18, 12))

# Plot 1: Scree plot
ax1 = plt.subplot2grid((2, 3), (0, 0))
ax1.bar(eigenvals['pc'][:10], eigenvals['variance_explained'][:10], 
        color='skyblue', edgecolor='black', alpha=0.8)
ax1.plot(eigenvals['pc'][:10], eigenvals['variance_explained'][:10], 
         'ro-', linewidth=2, markersize=6)
ax1.set_title('Scree Plot - Variance Explained', fontsize=14, fontweight='bold')
ax1.set_ylabel('Variance Explained (%)', fontsize=12)
ax1.set_xlabel('Principal Component', fontsize=12)
ax1.tick_params(axis='x', rotation=45)

# Add values on bars
for i, v in enumerate(eigenvals['variance_explained'][:10]):
    ax1.text(i, v + 0.5, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')

# Plot 2: PC1 vs PC2
ax2 = plt.subplot2grid((2, 3), (0, 1))
scatter = ax2.scatter(pca_df['PC1'], pca_df['PC2'], s=100, alpha=0.7, 
                     c=range(len(pca_df)), cmap='viridis')
ax2.set_title(f'PC1 vs PC2\n({eigenvals["variance_explained"].iloc[0]:.1f}% vs {eigenvals["variance_explained"].iloc[1]:.1f}%)', 
              fontsize=14, fontweight='bold')
ax2.set_xlabel(f'PC1 ({eigenvals["variance_explained"].iloc[0]:.1f}%)', fontsize=12)
ax2.set_ylabel(f'PC2 ({eigenvals["variance_explained"].iloc[1]:.1f}%)', fontsize=12)

# Add sample labels
for i, row in pca_df.iterrows():
    ax2.annotate(row['IID'], (row['PC1'], row['PC2']), 
                xytext=(5, 5), textcoords='offset points', 
                fontsize=8, alpha=0.8)

# Plot 3: PC1 vs PC3
ax3 = plt.subplot2grid((2, 3), (0, 2))
scatter = ax3.scatter(pca_df['PC1'], pca_df['PC3'], s=100, alpha=0.7, 
                     c=range(len(pca_df)), cmap='viridis')
ax3.set_title(f'PC1 vs PC3\n({eigenvals["variance_explained"].iloc[0]:.1f}% vs {eigenvals["variance_explained"].iloc[2]:.1f}%)', 
              fontsize=14, fontweight='bold')
ax3.set_xlabel(f'PC1 ({eigenvals["variance_explained"].iloc[0]:.1f}%)', fontsize=12)
ax3.set_ylabel(f'PC3 ({eigenvals["variance_explained"].iloc[2]:.1f}%)', fontsize=12)

for i, row in pca_df.iterrows():
    ax3.annotate(row['IID'], (row['PC1'], row['PC3']), 
                xytext=(5, 5), textcoords='offset points', 
                fontsize=8, alpha=0.8)

# Plot 4: PC2 vs PC3
ax4 = plt.subplot2grid((2, 3), (1, 0))
scatter = ax4.scatter(pca_df['PC2'], pca_df['PC3'], s=100, alpha=0.7, 
                     c=range(len(pca_df)), cmap='viridis')
ax4.set_title(f'PC2 vs PC3\n({eigenvals["variance_explained"].iloc[1]:.1f}% vs {eigenvals["variance_explained"].iloc[2]:.1f}%)', 
              fontsize=14, fontweight='bold')
ax4.set_xlabel(f'PC2 ({eigenvals["variance_explained"].iloc[1]:.1f}%)', fontsize=12)
ax4.set_ylabel(f'PC3 ({eigenvals["variance_explained"].iloc[2]:.1f}%)', fontsize=12)

for i, row in pca_df.iterrows():
    ax4.annotate(row['IID'], (row['PC2'], row['PC3']), 
                xytext=(5, 5), textcoords='offset points', 
                fontsize=8, alpha=0.8)

# Plot 5: Cumulative variance
ax5 = plt.subplot2grid((2, 3), (1, 1))
cumulative_variance = eigenvals['variance_explained'].cumsum()
ax5.plot(range(1, 11), cumulative_variance[:10], 'bo-', linewidth=2, markersize=6)
ax5.fill_between(range(1, 11), cumulative_variance[:10], alpha=0.3)
ax5.set_title('Cumulative Variance Explained', fontsize=14, fontweight='bold')
ax5.set_xlabel('Number of Principal Components', fontsize=12)
ax5.set_ylabel('Cumulative Variance (%)', fontsize=12)
ax5.grid(True, alpha=0.3)

# Add values on points
for i, v in enumerate(cumulative_variance[:10]):
    ax5.text(i+1, v + 2, f'{v:.1f}%', ha='center', va='bottom', fontweight='bold')

# Plot 6: PC loadings heatmap (first 5 PCs for top samples)
ax6 = plt.subplot2grid((2, 3), (1, 2))
# Use top 10 samples for clarity in heatmap
top_samples = pca_df.nlargest(10, 'PC1')['IID'].tolist() + pca_df.nsmallest(10, 'PC1')['IID'].tolist()
top_samples = list(dict.fromkeys(top_samples))  # Remove duplicates
heatmap_data = pca_df[pca_df['IID'].isin(top_samples)].set_index('IID')[['PC1', 'PC2', 'PC3', 'PC4', 'PC5']]

im = ax6.imshow(heatmap_data.T, cmap='coolwarm', aspect='auto', 
                vmin=-heatmap_data.abs().max().max(), 
                vmax=heatmap_data.abs().max().max())
ax6.set_title('PCA Loadings Heatmap\n(Top/Bottom Samples)', fontsize=14, fontweight='bold')
ax6.set_ylabel('Principal Components', fontsize=12)
ax6.set_xlabel('Samples', fontsize=12)
ax6.set_xticks(range(len(heatmap_data.index)))
ax6.set_xticklabels(heatmap_data.index, rotation=45, ha='right')
ax6.set_yticks(range(5))
ax6.set_yticklabels(['PC1', 'PC2', 'PC3', 'PC4', 'PC5'])

# Add colorbar
cbar = plt.colorbar(im, ax=ax6, shrink=0.8)
cbar.set_label('PCA Loading', rotation=270, labelpad=15)

plt.tight_layout()
plt.savefig('Diversity/PCA/Faba_PCA_Comprehensive.png', dpi=300, bbox_inches='tight')
plt.show()

print("✓ Comprehensive PCA analysis saved as 'Diversity/PCA/Faba_PCA_Comprehensive.png'")

# Save PCA summary
pca_summary = {
    'total_variance_top10': eigenvals['variance_explained'].sum(),
    'pc1_variance': eigenvals['variance_explained'].iloc[0],
    'pc2_variance': eigenvals['variance_explained'].iloc[1],
    'pc3_variance': eigenvals['variance_explained'].iloc[2],
    'n_components_90pct': (cumulative_variance >= 90).idxmax() + 1 if (cumulative_variance >= 90).any() else 10
}

with open('Diversity/PCA/PCA_summary_stats.txt', 'w') as f:
    f.write("PCA Analysis Summary Statistics\n")
    f.write("===============================\n")
    for key, value in pca_summary.items():
        f.write(f"{key}: {value:.2f}\n")
    
    f.write("\nVariance Explained by Each PC:\n")
    for i, row in eigenvals.iterrows():
        f.write(f"{row['pc']}: {row['variance_explained']:.2f}%\n")

print("✓ PCA summary statistics saved as 'Diversity/PCA/PCA_summary_stats.txt'")

# Print population structure insights
print(f"\n🔍 POPULATION STRUCTURE INSIGHTS:")
print(f"PC1 explains {eigenvals['variance_explained'].iloc[0]:.1f}% of variance")
print(f"PC2 explains {eigenvals['variance_explained'].iloc[1]:.1f}% of variance")
print(f"PC3 explains {eigenvals['variance_explained'].iloc[2]:.1f}% of variance")
print(f"Top 3 PCs explain {cumulative_variance.iloc[2]:.1f}% of total variance")
print(f"Components needed for 90% variance: {pca_summary['n_components_90pct']}")

# Identify potential subgroups based on PC1
pc1_threshold = pca_df['PC1'].median()
group1 = pca_df[pca_df['PC1'] > pc1_threshold]['IID'].tolist()
group2 = pca_df[pca_df['PC1'] <= pc1_threshold]['IID'].tolist()

print(f"\n📊 POTENTIAL SUBGROUPS (based on PC1 median split):")
print(f"Group 1 (High PC1, {len(group1)} samples): {', '.join(group1)}")
print(f"Group 2 (Low PC1, {len(group2)} samples): {', '.join(group2)}")
