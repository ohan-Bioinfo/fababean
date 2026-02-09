# PhylogeneticTree/scripts/generate_distance_matrix.py
import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from Bio.Phylo.TreeConstruction import DistanceMatrix
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
import subprocess

def calculate_distance_matrices():
    """Calculate various distance matrices for phylogenetic analysis"""
    
    # Read numerical genotypes
    raw = pd.read_csv('data/faba_fingerprint.raw', sep='\s+')
    sample_ids = raw['IID'].tolist()
    
    # Extract SNP columns
    snp_cols = [col for col in raw.columns if col.startswith('SNP')]
    genotype_matrix = raw[snp_cols].values
    
    print(f"Calculating distance matrices for {len(sample_ids)} samples...")
    
    # 1. Euclidean distance
    euclidean_dist = pdist(genotype_matrix, metric='euclidean')
    euclidean_matrix = squareform(euclidean_dist)
    pd.DataFrame(euclidean_matrix, index=sample_ids, columns=sample_ids)\
      .to_csv('output/euclidean_distance_matrix.csv')
    
    # 2. Manhattan distance
    manhattan_dist = pdist(genotype_matrix, metric='cityblock')
    manhattan_matrix = squareform(manhattan_dist)
    pd.DataFrame(manhattan_matrix, index=sample_ids, columns=sample_ids)\
      .to_csv('output/manhattan_distance_matrix.csv')
    
    # 3. Hamming distance (for genetic data)
    hamming_dist = pdist(genotype_matrix, metric='hamming')
    hamming_matrix = squareform(hamming_dist)
    pd.DataFrame(hamming_matrix, index=sample_ids, columns=sample_ids)\
      .to_csv('output/hamming_distance_matrix.csv')
    
    # 4. Create PHYLIP format distance matrix
    create_phylip_distance_matrix(sample_ids, euclidean_matrix, 'output/euclidean_dist.phy')
    create_phylip_distance_matrix(sample_ids, hamming_matrix, 'output/hamming_dist.phy')
    
    print("Distance matrices calculated and saved")
    return sample_ids, euclidean_matrix, hamming_matrix

def create_phylip_distance_matrix(sample_ids, distance_matrix, filename):
    """Create PHYLIP format distance matrix"""
    n_samples = len(sample_ids)
    
    with open(filename, 'w') as f:
        f.write(f' {n_samples}\n')
        for i, sample_id in enumerate(sample_ids):
            # Format sample ID (max 10 characters)
            formatted_id = sample_id[:10].ljust(10)
            distances = ' '.join([f'{dist:.6f}' for dist in distance_matrix[i]])
            f.write(f'{formatted_id} {distances}\n')
    
    print(f"PHYLIP distance matrix saved: {filename}")

if __name__ == "__main__":
    sample_ids, euclidean_matrix, hamming_matrix = calculate_distance_matrices()