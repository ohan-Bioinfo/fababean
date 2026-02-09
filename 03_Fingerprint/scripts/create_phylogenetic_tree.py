import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import subprocess
import os
from matplotlib import rcParams

# Set style
plt.style.use('default')
rcParams['font.family'] = 'DejaVu Sans'

def create_phylogenetic_tree():
    print("=== CONSTRUCTING PHYLOGENETIC TREE USING 150 SNPS ===")
    
    # Check if we have the necessary files
    required_files = [
        'data/top_150_snps_list.txt',
        'data/Faba_high_quality.bed',
        'data/Faba_high_quality.bim', 
        'data/Faba_high_quality.fam'
    ]
    
    for file in required_files:
        if not os.path.exists(file):
            print(f"❌ Missing file: {file}")
            return
    
    print("✓ All required files found")
    
    # Step 1: Extract the 150 SNPs using PLINK
    print("\n1. Extracting 150 SNPs using PLINK...")
    
    plink_cmd = [
        "plink",
        "--bfile", "data/Faba_high_quality",
        "--extract", "data/top_150_snps_list.txt",
        "--make-bed",
        "--out", "output/faba_150_snps_tree",
        "--allow-extra-chr"
    ]
    
    result = subprocess.run(plink_cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print("✓ Successfully extracted 150 SNPs")
    else:
        print(f"❌ PLINK error: {result.stderr}")
        return
    
    # Step 2: Convert to VCF for easier processing
    print("\n2. Converting to VCF format...")
    
    plink_cmd = [
        "plink",
        "--bfile", "output/faba_150_snps_tree",
        "--recode", "vcf",
        "--out", "output/faba_150_snps_tree",
        "--allow-extra-chr"
    ]
    
    result = subprocess.run(plink_cmd, capture_output=True, text=True)
    if result.returncode == 0:
        print("✓ Successfully converted to VCF")
    else:
        print(f"⚠️ VCF conversion warning: {result.stderr}")
    
    # Step 3: Read and process the genotype data
    print("\n3. Processing genotype data...")
    
    # Read the BIM file to get SNP information
    bim_data = pd.read_csv('output/faba_150_snps_tree.bim', 
                          sep='\t', 
                          header=None,
                          names=['chr', 'snp_id', 'cM', 'position', 'allele1', 'allele2'])
    
    # Read the FAM file to get sample information
    fam_data = pd.read_csv('output/faba_150_snps_tree.fam', 
                          sep='\t', 
                          header=None,
                          names=['family', 'sample', 'father', 'mother', 'sex', 'phenotype'])
    
    samples = fam_data['sample'].tolist()
    print(f"✓ Loaded {len(samples)} samples: {samples}")
    print(f"✓ Using {len(bim_data)} SNPs")
    
    # Step 4: Read the RAW file to get genotype matrix
    print("\n4. Reading genotype data...")
    
    plink_cmd = [
        "plink",
        "--bfile", "output/faba_150_snps_tree",
        "--recode", "A",
        "--out", "output/faba_150_snps_tree",
        "--allow-extra-chr"
    ]
    
    result = subprocess.run(plink_cmd, capture_output=True, text=True)
    
    if os.path.exists('output/faba_150_snps_tree.raw'):
        raw_data = pd.read_csv('output/faba_150_snps_tree.raw', sep='\s+')
        print(f"✓ Loaded genotype data for {len(raw_data)} samples")
        
        # Extract genotype matrix (remove first 6 columns which are metadata)
        genotype_columns = [col for col in raw_data.columns if col not in ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        genotype_matrix = raw_data[genotype_columns].values
        
        print(f"✓ Genotype matrix shape: {genotype_matrix.shape}")
        
    else:
        print("❌ Could not read genotype data")
        return
    
    # Step 5: Calculate genetic distance matrix
    print("\n5. Calculating genetic distance matrix...")
    
    # Handle missing data (coded as -9 in PLINK)
    genotype_matrix[genotype_matrix == -9] = np.nan
    
    # Calculate pairwise genetic distances (Euclidean distance)
    # Transpose to get samples as rows, SNPs as columns
    distance_matrix = pdist(genotype_matrix, metric='euclidean')
    distance_square = squareform(distance_matrix)
    
    print(f"✓ Distance matrix calculated")
    print(f"  Distance range: {np.nanmin(distance_square):.2f} - {np.nanmax(distance_square):.2f}")
    
    # Step 6: Build phylogenetic tree using hierarchical clustering
    print("\n6. Building phylogenetic tree...")
    
    # Perform hierarchical clustering
    linkage_matrix = linkage(distance_matrix, method='average')
    
    # Step 7: Create the phylogenetic tree visualization
    print("\n7. Creating phylogenetic tree visualization...")
    
    # Create figure with multiple views
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(20, 16))
    
    # Plot 1: Standard dendrogram
    dendrogram(linkage_matrix, 
               labels=samples,
               orientation='top',
               leaf_rotation=45,
               ax=ax1)
    ax1.set_title('Phylogenetic Tree - 150 SNPs\n(Hierarchical Clustering)', 
                  fontsize=14, fontweight='bold', pad=20)
    ax1.set_xlabel('Accession', fontsize=12, labelpad=10)
    ax1.set_ylabel('Genetic Distance', fontsize=12, labelpad=10)
    
    # Plot 2: Left-oriented dendrogram
    dendrogram(linkage_matrix, 
               labels=samples,
               orientation='left',
               leaf_rotation=0,
               ax=ax2)
    ax2.set_title('Phylogenetic Tree\n(Left Orientation)', 
                  fontsize=14, fontweight='bold', pad=20)
    ax2.set_xlabel('Genetic Distance', fontsize=12, labelpad=10)
    
    # Plot 3: Heatmap of genetic distances
    sns.heatmap(distance_square, 
                xticklabels=samples,
                yticklabels=samples,
                cmap='viridis',
                square=True,
                ax=ax3)
    ax3.set_title('Genetic Distance Matrix\n(150 SNPs)', 
                  fontsize=14, fontweight='bold', pad=20)
    ax3.set_xlabel('Accession', fontsize=12, labelpad=10)
    ax3.set_ylabel('Accession', fontsize=12, labelpad=10)
    
    # Plot 4: Clustered heatmap
    sns.clustermap(pd.DataFrame(distance_square, index=samples, columns=samples),
                  cmap='viridis',
                  row_linkage=linkage_matrix,
                  col_linkage=linkage_matrix,
                  xticklabels=True,
                  yticklabels=True,
                  ax=ax4)
    ax4.set_title('Clustered Distance Matrix', 
                  fontsize=14, fontweight='bold', pad=20)
    
    plt.tight_layout()
    plt.savefig('plots/phylogenetic_tree_150_snps.png', dpi=300, bbox_inches='tight')
    plt.savefig('plots/phylogenetic_tree_150_snps.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Step 8: Create a clean standalone tree
    print("\n8. Creating clean standalone tree...")
    
    plt.figure(figsize=(12, 8))
    dendrogram(linkage_matrix, 
               labels=samples,
               orientation='right',
               leaf_rotation=0,
               leaf_font_size=10,
               color_threshold=0.7 * max(linkage_matrix[:, 2]))
    
    plt.title('Phylogenetic Tree of Faba Bean Accessions\n(150 Informative SNPs)', 
              fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Genetic Distance', fontsize=12, labelpad=10)
    
    plt.tight_layout()
    plt.savefig('plots/phylogenetic_tree_clean.png', dpi=300, bbox_inches='tight')
    plt.savefig('plots/phylogenetic_tree_clean.pdf', dpi=300, bbox_inches='tight')
    plt.show()
    
    # Step 9: Save cluster assignments
    print("\n9. Saving cluster assignments...")
    
    # Cut the tree to get clusters (you can adjust the threshold)
    clusters = fcluster(linkage_matrix, t=3, criterion='maxclust')
    
    cluster_assignments = pd.DataFrame({
        'Sample': samples,
        'Cluster': clusters,
        'Distance_To_Centroid': [np.mean(distance_square[i, clusters == clusters[i]]) 
                                for i in range(len(samples))]
    })
    
    cluster_assignments.to_csv('output/phylogenetic_clusters.csv', index=False)
    
    # Step 10: Print summary statistics
    print("\n📊 PHYLOGENETIC TREE SUMMARY:")
    print(f"✓ Samples: {len(samples)}")
    print(f"✓ SNPs used: {len(bim_data)}")
    print(f"✓ Genetic distance range: {np.min(distance_square):.2f} - {np.max(distance_square):.2f}")
    print(f"✓ Number of clusters identified: {len(np.unique(clusters))}")
    
    print(f"\n🔍 CLUSTER ASSIGNMENTS:")
    for cluster_num in np.unique(clusters):
        cluster_samples = cluster_assignments[cluster_assignments['Cluster'] == cluster_num]['Sample'].tolist()
        print(f"  Cluster {cluster_num}: {len(cluster_samples)} samples")
        print(f"    {cluster_samples}")
    
    # Find most similar and most distant pairs
    print(f"\n🔗 RELATIONSHIP ANALYSIS:")
    
    # Most similar pairs (excluding self-comparisons)
    np.fill_diagonal(distance_square, np.inf)  # Ignore diagonal
    min_idx = np.unravel_index(np.argmin(distance_square), distance_square.shape)
    min_distance = distance_square[min_idx]
    print(f"  Most similar pair: {samples[min_idx[0]]} - {samples[min_idx[1]]} (distance: {min_distance:.3f})")
    
    # Most distant pairs
    max_idx = np.unravel_index(np.argmax(distance_square), distance_square.shape)
    max_distance = distance_square[max_idx]
    print(f"  Most distant pair: {samples[max_idx[0]]} - {samples[max_idx[1]]} (distance: {max_distance:.3f})")
    
    # Save distance matrix
    distance_df = pd.DataFrame(distance_square, index=samples, columns=samples)
    distance_df.to_csv('output/genetic_distance_matrix.csv')
    
    print(f"\n💾 FILES SAVED:")
    print(f"  - plots/phylogenetic_tree_150_snps.png/pdf")
    print(f"  - plots/phylogenetic_tree_clean.png/pdf") 
    print(f"  - output/phylogenetic_clusters.csv")
    print(f"  - output/genetic_distance_matrix.csv")
    
    return cluster_assignments, distance_df

def create_alternative_tree():
    """Alternative method using different distance metrics"""
    print("\n=== CREATING ALTERNATIVE TREE (Hamming Distance) ===")
    
    # Read the RAW file again for alternative processing
    if os.path.exists('output/faba_150_snps_tree.raw'):
        raw_data = pd.read_csv('output/faba_150_snps_tree.raw', sep='\s+')
        samples = raw_data['IID'].tolist()
        
        # Extract genotype matrix
        genotype_columns = [col for col in raw_data.columns if col not in ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']]
        genotype_matrix = raw_data[genotype_columns].values
        
        # Convert to binary (0/1) for Hamming distance
        # Handle missing data and convert to 0/1/2 format
        binary_matrix = np.where(genotype_matrix == -9, np.nan, genotype_matrix)
        
        # Calculate Hamming distance (proportion of differing SNPs)
        from sklearn.metrics.pairwise import pairwise_distances
        hamming_distances = pairwise_distances(binary_matrix, metric='hamming')
        
        # Build tree
        linkage_hamming = linkage(squareform(hamming_distances), method='ward')
        
        # Plot
        plt.figure(figsize=(12, 8))
        dendrogram(linkage_hamming, 
                  labels=samples,
                  orientation='right',
                  leaf_rotation=0)
        
        plt.title('Phylogenetic Tree - Hamming Distance\n(150 SNPs)', 
                  fontsize=16, fontweight='bold', pad=20)
        plt.xlabel('Hamming Distance', fontsize=12, labelpad=10)
        
        plt.tight_layout()
        plt.savefig('plots/phylogenetic_tree_hamming.png', dpi=300, bbox_inches='tight')
        plt.show()
        
        print("✓ Alternative tree (Hamming distance) created")

if __name__ == "__main__":
    # Create main phylogenetic tree
    clusters, distances = create_phylogenetic_tree()
    
    # Create alternative tree with different distance metric
    create_alternative_tree()
