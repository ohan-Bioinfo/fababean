import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
import pandas as pd

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

def create_simple_side_by_side():
    print("=== SIDE-BY-SIDE KINSHIP & IBS HEATMAPS ===")
    
    # Read sample IDs
    with open('Diversity/Kinship/Faba_Kinship.king.id', 'r') as f:
        lines = f.readlines()
        samples = [line.strip().split()[1] for line in lines[1:]]

    print(f"Number of samples: {len(samples)}")

    # Read KINSHIP data
    kinship_lines = []
    with open('Diversity/Kinship/Faba_Kinship.king', 'r') as f:
        for line in f:
            if line.strip():
                kinship_lines.append([float(x) for x in line.strip().split()])

    n = len(samples)
    kinship_matrix = np.zeros((n, n))
    
    for i in range(len(kinship_lines)):
        for j in range(len(kinship_lines[i])):
            kinship_matrix[i+1, j] = kinship_lines[i][j]
            kinship_matrix[j, i+1] = kinship_lines[i][j]
    
    np.fill_diagonal(kinship_matrix, 1.0)

    # Read IBS data
    ibs_matrix = np.loadtxt('Diversity/IBS/Faba_IBS.mibs')

    # Create simple side-by-side figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # Kinship heatmap
    im1 = ax1.imshow(kinship_matrix, cmap='RdYlBu_r', aspect='auto', 
                    vmin=-0.5, vmax=0.5)
    ax1.set_xticks(range(n))
    ax1.set_yticks(range(n))
    ax1.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax1.set_yticklabels(samples, fontsize=10)
    ax1.set_xlabel('Accession', fontsize=12, labelpad=10)
    ax1.set_ylabel('Accession', fontsize=12, labelpad=10)
    ax1.set_title('KINSHIP Matrix\nFaba Bean Germplasm', 
                  fontsize=14, fontweight='bold', pad=15)
    cbar1 = plt.colorbar(im1, ax=ax1, shrink=0.8)
    cbar1.set_label('Kinship Coefficient', rotation=270, labelpad=20)

    # IBS heatmap
    im2 = ax2.imshow(ibs_matrix, cmap='RdYlBu_r', aspect='auto', 
                    vmin=0.75, vmax=1.0)
    ax2.set_xticks(range(n))
    ax2.set_yticks(range(n))
    ax2.set_xticklabels(samples, rotation=45, ha='right', fontsize=10)
    ax2.set_yticklabels(samples, fontsize=10)
    ax2.set_xlabel('Accession', fontsize=12, labelpad=10)
    ax2.set_ylabel('Accession', fontsize=12, labelpad=10)
    ax2.set_title('IBS Similarity Matrix\nFaba Bean Germplasm', 
                  fontsize=14, fontweight='bold', pad=15)
    cbar2 = plt.colorbar(im2, ax=ax2, shrink=0.8)
    cbar2.set_label('IBS Proportion', rotation=270, labelpad=20)

    # Add grid to both heatmaps
    for ax in [ax1, ax2]:
        ax.set_xticks(np.arange(-0.5, n, 1), minor=True)
        ax.set_yticks(np.arange(-0.5, n, 1), minor=True)
        ax.grid(which="minor", color="white", linestyle='-', linewidth=0.5)
        ax.tick_params(which="minor", size=0)

    plt.tight_layout()
    plt.savefig('Diversity/SideBySide_Kinship_IBS.png', dpi=300, bbox_inches='tight')
    plt.savefig('Diversity/SideBySide_Kinship_IBS.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    print("✓ Side-by-side heatmaps saved as:")
    print("  - Diversity/SideBySide_Kinship_IBS.png")
    print("  - Diversity/SideBySide_Kinship_IBS.pdf")

# Run the function
create_simple_side_by_side()