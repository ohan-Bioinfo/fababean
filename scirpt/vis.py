import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import rcParams
from matplotlib.patches import FancyBboxPatch
import os

# Set style for better plots
plt.style.use('seaborn-v0_8-whitegrid')
rcParams['font.family'] = 'DejaVu Sans'
rcParams['font.size'] = 10
rcParams['figure.figsize'] = (14, 10)

def load_pca_data(eigenvec_file, eigenval_file):
    """Load PCA results from eigenvec and eigenval files"""
    # Load eigenvectors
    eigenvec = pd.read_csv(eigenvec_file, sep='\s+', header=None)
    
    # Set column names - first two columns are sample IDs, rest are PCs
    pc_columns = [f'PC{i+1}' for i in range(len(eigenvec.columns)-2)]
    eigenvec.columns = ['FID', 'IID'] + pc_columns
    
    # Load eigenvalues
    eigenval = pd.read_csv(eigenval_file, sep='\s+', header=None, names=['Eigenvalue'])
    
    return eigenvec, eigenval

def calculate_variance_explained(eigenval):
    """Calculate variance explained by each principal component"""
    total_variance = eigenval['Eigenvalue'].sum()
    eigenval['Variance_Explained'] = (eigenval['Eigenvalue'] / total_variance) * 100
    eigenval['Cumulative_Variance'] = eigenval['Variance_Explained'].cumsum()
    eigenval['PC'] = [f'PC{i+1}' for i in range(len(eigenval))]
    return eigenval

def create_pca_summary(eigenvec, eigenval):
    """Create comprehensive PCA summary"""
    summary = {
        'Total_Samples': len(eigenvec),
        'Total_PCs': len(eigenval),
        'Total_Variance': eigenval['Eigenvalue'].sum(),
        'PCs_for_80%_variance': len(eigenval[eigenval['Cumulative_Variance'] <= 80]),
        'PCs_for_90%_variance': len(eigenval[eigenval['Cumulative_Variance'] <= 90]),
        'PC1_Variance': eigenval.iloc[0]['Variance_Explained'],
        'PC2_Variance': eigenval.iloc[1]['Variance_Explained'],
        'PC3_Variance': eigenval.iloc[2]['Variance_Explained']
    }
    return summary

def add_sample_labels(ax, x, y, labels, line_alpha=0.3, text_alpha=0.8):
    """Add sample labels with light connecting lines to avoid overlapping"""
    for i, (xi, yi, label) in enumerate(zip(x, y, labels)):
        # Add a light gray line from point to text
        ax.plot([xi, xi], [yi, yi], 'gray', alpha=line_alpha, linewidth=0.5)
        
        # Add text label with slight offset
        ax.annotate(label, (xi, yi), 
                   xytext=(5, 5), textcoords='offset points',
                   fontsize=8, alpha=text_alpha,
                   bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                           alpha=0.7, edgecolor='none'),
                   arrowprops=dict(arrowstyle='-', color='gray', 
                                 alpha=line_alpha, linewidth=0.5))

def plot_pca_comprehensive(eigenvec, eigenval, output_dir):
    """Create comprehensive PCA visualization with sample labels"""
    
    # Create subplots
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Plot 1: PC1 vs PC2 with labels
    ax1 = axes[0, 0]
    scatter1 = ax1.scatter(eigenvec['PC1'], eigenvec['PC2'], alpha=0.8, 
                          c=np.arange(len(eigenvec)), cmap='viridis', s=80,
                          edgecolor='black', linewidth=0.5)
    ax1.set_xlabel(f'PC1 ({eigenval.iloc[0]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax1.set_ylabel(f'PC2 ({eigenval.iloc[1]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax1.set_title('PCA: PC1 vs PC2 with Sample Labels', fontsize=14, fontweight='bold')
    
    # Add sample labels
    add_sample_labels(ax1, eigenvec['PC1'], eigenvec['PC2'], eigenvec['IID'])
    plt.colorbar(scatter1, ax=ax1, label='Sample Index')
    
    # Plot 2: PC1 vs PC3 with labels
    ax2 = axes[0, 1]
    scatter2 = ax2.scatter(eigenvec['PC1'], eigenvec['PC3'], alpha=0.8,
                          c=np.arange(len(eigenvec)), cmap='plasma', s=80,
                          edgecolor='black', linewidth=0.5)
    ax2.set_xlabel(f'PC1 ({eigenval.iloc[0]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax2.set_ylabel(f'PC3 ({eigenval.iloc[2]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax2.set_title('PCA: PC1 vs PC3 with Sample Labels', fontsize=14, fontweight='bold')
    
    # Add sample labels
    add_sample_labels(ax2, eigenvec['PC1'], eigenvec['PC3'], eigenvec['IID'])
    plt.colorbar(scatter2, ax=ax2, label='Sample Index')
    
    # Plot 3: Scree plot
    ax3 = axes[1, 0]
    pcs_to_plot = min(10, len(eigenval))
    bars = ax3.bar(range(1, pcs_to_plot+1), 
                   eigenval['Variance_Explained'].head(pcs_to_plot),
                   color='skyblue', alpha=0.7, edgecolor='navy', linewidth=0.5)
    ax3.set_xlabel('Principal Component', fontsize=12)
    ax3.set_ylabel('Variance Explained (%)', fontsize=12)
    ax3.set_title('Scree Plot - Variance per Principal Component', fontsize=14, fontweight='bold')
    ax3.set_xticks(range(1, pcs_to_plot+1))
    
    # Add variance percentage on bars
    for bar, variance in zip(bars, eigenval['Variance_Explained'].head(pcs_to_plot)):
        height = bar.get_height()
        ax3.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                f'{variance:.1f}%', ha='center', va='bottom', fontsize=9, fontweight='bold')
    
    # Plot 4: Cumulative variance
    ax4 = axes[1, 1]
    ax4.plot(range(1, pcs_to_plot+1), 
             eigenval['Cumulative_Variance'].head(pcs_to_plot), 
             'bo-', linewidth=2, markersize=6, markerfacecolor='blue')
    ax4.axhline(y=80, color='red', linestyle='--', alpha=0.7, linewidth=1.5, label='80% variance')
    ax4.axhline(y=90, color='green', linestyle='--', alpha=0.7, linewidth=1.5, label='90% variance')
    ax4.set_xlabel('Number of Principal Components', fontsize=12)
    ax4.set_ylabel('Cumulative Variance Explained (%)', fontsize=12)
    ax4.set_title('Cumulative Variance Explained', fontsize=14, fontweight='bold')
    ax4.legend(fontsize=10)
    ax4.set_xticks(range(1, pcs_to_plot+1))
    ax4.grid(True, alpha=0.3)
    
    # Add value annotations for cumulative variance
    for i, (x_val, y_val) in enumerate(zip(range(1, pcs_to_plot+1), 
                                          eigenval['Cumulative_Variance'].head(pcs_to_plot))):
        ax4.annotate(f'{y_val:.1f}%', (x_val, y_val), 
                    xytext=(0, 10), textcoords='offset points',
                    ha='center', va='bottom', fontsize=8,
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='yellow', alpha=0.3))
    
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Faba_PCA_Comprehensive_Labeled.png'), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()
    
    return fig

def plot_pca_2d_labeled(eigenvec, eigenval, output_dir):
    """Create a detailed 2D PCA plot with sample labels"""
    plt.figure(figsize=(12, 9))
    
    # Create scatter plot
    scatter = plt.scatter(eigenvec['PC1'], eigenvec['PC2'], 
                         alpha=0.8, s=100, c=np.arange(len(eigenvec)), 
                         cmap='viridis', edgecolor='black', linewidth=0.8)
    
    plt.xlabel(f'Principal Component 1 ({eigenval.iloc[0]["Variance_Explained"]:.2f}%)', fontsize=14)
    plt.ylabel(f'Principal Component 2 ({eigenval.iloc[1]["Variance_Explained"]:.2f}%)', fontsize=14)
    plt.title('Faba Bean PCA - PC1 vs PC2 with Sample Labels', fontsize=16, fontweight='bold')
    
    # Add sample labels with light connecting lines
    add_sample_labels(plt.gca(), eigenvec['PC1'], eigenvec['PC2'], eigenvec['IID'])
    
    # Add colorbar
    cbar = plt.colorbar(scatter)
    cbar.set_label('Sample Index', fontsize=12)
    
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Faba_PCA_2D_Labeled.png'), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

def plot_pca_3d_labeled(eigenvec, eigenval, output_dir):
    """Create 3D PCA plot with sample labels"""
    from mpl_toolkits.mplot3d import Axes3D
    
    fig = plt.figure(figsize=(12, 9))
    ax = fig.add_subplot(111, projection='3d')
    
    scatter = ax.scatter(eigenvec['PC1'], eigenvec['PC2'], eigenvec['PC3'],
                        c=np.arange(len(eigenvec)), cmap='viridis', 
                        s=80, alpha=0.8, edgecolor='black', linewidth=0.5)
    
    ax.set_xlabel(f'PC1 ({eigenval.iloc[0]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax.set_ylabel(f'PC2 ({eigenval.iloc[1]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax.set_zlabel(f'PC3 ({eigenval.iloc[2]["Variance_Explained"]:.2f}%)', fontsize=12)
    ax.set_title('3D PCA - Faba Bean Samples', fontsize=16, fontweight='bold')
    
    # Add sample labels in 3D
    for i, (x, y, z, label) in enumerate(zip(eigenvec['PC1'], eigenvec['PC2'], 
                                            eigenvec['PC3'], eigenvec['IID'])):
        ax.text(x, y, z, label, fontsize=8, 
               bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    plt.colorbar(scatter, ax=ax, label='Sample Index', shrink=0.7)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, 'Faba_PCA_3D_Labeled.png'), 
                dpi=300, bbox_inches='tight', facecolor='white')
    plt.show()

def save_summary_files(eigenval, summary, eigenvec, output_dir):
    """Save summary statistics to files"""
    
    # Save detailed variance information
    variance_summary = eigenval[['PC', 'Eigenvalue', 'Variance_Explained', 'Cumulative_Variance']]
    variance_summary.to_csv(os.path.join(output_dir, 'PCA_variance_summary.csv'), index=False)
    
    # Save PC coordinates with sample IDs
    pc_coordinates = eigenvec.copy()
    pc_coordinates.to_csv(os.path.join(output_dir, 'PCA_sample_coordinates.csv'), index=False)
    
    # Save overall summary
    with open(os.path.join(output_dir, 'PCA_summary_statistics.txt'), 'w') as f:
        f.write("FABA BEAN PCA ANALYSIS SUMMARY\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"Total number of samples: {summary['Total_Samples']}\n")
        f.write(f"Total number of PCs: {summary['Total_PCs']}\n")
        f.write(f"Total variance: {summary['Total_Variance']:.2f}\n\n")
        
        f.write("Variance Explained by Top PCs:\n")
        f.write(f"  PC1: {summary['PC1_Variance']:.2f}%\n")
        f.write(f"  PC2: {summary['PC2_Variance']:.2f}%\n")
        f.write(f"  PC3: {summary['PC3_Variance']:.2f}%\n\n")
        
        f.write("Cumulative Variance Thresholds:\n")
        f.write(f"  PCs needed for 80% variance: {summary['PCs_for_80%_variance']}\n")
        f.write(f"  PCs needed for 90% variance: {summary['PCs_for_90%_variance']}\n\n")
        
        f.write("Sample IDs:\n")
        for sample_id in eigenvec['IID']:
            f.write(f"  {sample_id}\n")
        f.write("\n")
        
        f.write("Top 10 Principal Components:\n")
        for i in range(min(10, len(eigenval))):
            f.write(f"  {eigenval['PC'].iloc[i]}: {eigenval['Variance_Explained'].iloc[i]:.2f}% "
                   f"(Cumulative: {eigenval['Cumulative_Variance'].iloc[i]:.2f}%)\n")

def main():
    """Main function to run PCA analysis"""
    
    # Set paths
    pca_dir = "Diversity/PCA/"
    eigenvec_file = os.path.join(pca_dir, "Faba_PCA.eigenvec")
    eigenval_file = os.path.join(pca_dir, "Faba_PCA.eigenval")
    
    # Check if files exist
    if not os.path.exists(eigenvec_file):
        print(f"Error: {eigenvec_file} not found!")
        return
    if not os.path.exists(eigenval_file):
        print(f"Error: {eigenval_file} not found!")
        return
    
    print("Loading PCA data...")
    eigenvec, eigenval = load_pca_data(eigenvec_file, eigenval_file)
    
    print("Calculating variance explained...")
    eigenval = calculate_variance_explained(eigenval)
    
    print("Creating PCA summary...")
    summary = create_pca_summary(eigenvec, eigenval)
    
    # Print quick summary to console
    print("\n" + "="*60)
    print("FABA BEAN PCA - QUICK SUMMARY")
    print("="*60)
    print(f"Samples: {summary['Total_Samples']}")
    print(f"Sample IDs: {', '.join(eigenvec['IID'].tolist())}")
    print(f"PCs: {summary['Total_PCs']}")
    print(f"PC1 variance: {summary['PC1_Variance']:.2f}%")
    print(f"PC2 variance: {summary['PC2_Variance']:.2f}%")
    print(f"PC3 variance: {summary['PC3_Variance']:.2f}%")
    print(f"PCs for 80% variance: {summary['PCs_for_80%_variance']}")
    print(f"PCs for 90% variance: {summary['PCs_for_90%_variance']}")
    
    # Create visualizations with labels
    print("\nCreating labeled visualizations...")
    plot_pca_comprehensive(eigenvec, eigenval, pca_dir)
    plot_pca_2d_labeled(eigenvec, eigenval, pca_dir)
    
    # Create 3D plot if you have at least 3 PCs
    if 'PC3' in eigenvec.columns:
        print("Creating 3D visualization...")
        plot_pca_3d_labeled(eigenvec, eigenval, pca_dir)
    
    # Save summary files
    print("Saving summary files...")
    save_summary_files(eigenval, summary, eigenvec, pca_dir)
    
    print(f"\nAnalysis complete! Check the {pca_dir} directory for:")
    print("  - Faba_PCA_Comprehensive_Labeled.png (4-panel plot with labels)")
    print("  - Faba_PCA_2D_Labeled.png (detailed 2D plot with labels)")
    print("  - Faba_PCA_3D_Labeled.png (3D PCA plot)")
    print("  - PCA_variance_summary.csv (detailed variance table)")
    print("  - PCA_sample_coordinates.csv (sample PC coordinates)")
    print("  - PCA_summary_statistics.txt (comprehensive summary)")

if __name__ == "__main__":
    main()