import matplotlib.pyplot as plt
import seaborn as sns
from data_loader import SNPDataLoader

def plot_missingness_distribution(data):
    """Plot sample and SNP missingness distributions"""
    
    if 'sample_missingness' not in data or 'snp_missingness' not in data:
        print("Missingness data not available")
        return
    
    fig, axes = plt.subplots(1, 2, figsize=(16, 6))
    fig.suptitle('Missingness Distribution', fontsize=16, fontweight='bold')
    
    # Sample missingness
    axes[0].hist(data['sample_missingness']['F_MISS'], bins=50, 
                alpha=0.7, color='coral', edgecolor='black')
    axes[0].axvline(0.05, color='red', linestyle='--', linewidth=2, 
                   label='5% threshold (typical)')
    axes[0].set_title('A. Sample Missingness', fontweight='bold')
    axes[0].set_xlabel('Missing Rate (F_MISS)')
    axes[0].set_ylabel('Number of Samples')
    axes[0].legend()
    axes[0].grid(alpha=0.3)
    
    # SNP missingness
    axes[1].hist(data['snp_missingness']['F_MISS'], bins=50, 
                alpha=0.7, color='lightseagreen', edgecolor='black')
    axes[1].axvline(0.05, color='red', linestyle='--', linewidth=2, 
                   label='5% threshold (typical)')
    axes[1].set_title('B. SNP Missingness', fontweight='bold')
    axes[1].set_xlabel('Missing Rate (F_MISS)')
    axes[1].set_ylabel('Number of SNPs')
    axes[1].legend()
    axes[1].grid(alpha=0.3)
    
    plt.tight_layout()
    plt.savefig('Faba_Bean_Missingness_Distribution.png', dpi=300, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    plot_missingness_distribution(data)
