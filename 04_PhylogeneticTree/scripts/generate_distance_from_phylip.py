# scripts/generate_distance_from_phylip.py
import pandas as pd
import numpy as np
import os

def read_phylip_file(phy_file):
    """Read PHYLIP format file and return sequence data"""
    print(f"Reading PHYLIP file: {phy_file}")
    
    if not os.path.exists(phy_file):
        print(f"ERROR: File not found: {phy_file}")
        return None, None, None
    
    with open(phy_file, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]
    
    if not lines:
        print("ERROR: PHYLIP file is empty")
        return None, None, None
    
    # First line: number of sequences and sequence length
    try:
        header_parts = lines[0].split()
        n_seqs = int(header_parts[0])
        seq_len = int(header_parts[1])
        print(f"PHYLIP header: {n_seqs} sequences, {seq_len} sites")
    except:
        print("ERROR: Invalid PHYLIP header format")
        return None, None, None
    
    sequences = {}
    for line in lines[1:]:
        if line:
            parts = line.split()
            if len(parts) >= 2:
                seq_id = parts[0]
                seq_data = ''.join(parts[1:])
                sequences[seq_id] = seq_data
                print(f"  Loaded sequence: {seq_id} (length: {len(seq_data)})")
    
    print(f"Successfully loaded {len(sequences)} sequences")
    return sequences, n_seqs, seq_len

def calculate_simple_distance_matrix(sequences):
    """Calculate simple Hamming distance matrix"""
    sample_ids = list(sequences.keys())
    n_samples = len(sample_ids)
    
    print("Calculating distance matrix...")
    
    # Initialize distance matrix
    distance_matrix = np.zeros((n_samples, n_samples))
    
    for i in range(n_samples):
        seq1 = sequences[sample_ids[i]]
        for j in range(i, n_samples):
            if i == j:
                distance_matrix[i, j] = 0
            else:
                seq2 = sequences[sample_ids[j]]
                # Calculate proportion of differing sites
                matches = sum(1 for a, b in zip(seq1, seq2) if a == b and a != '?' and a != '-')
                total_sites = sum(1 for a, b in zip(seq1, seq2) if a != '?' and b != '?' and a != '-' and b != '-')
                
                if total_sites > 0:
                    distance_matrix[i, j] = 1 - (matches / total_sites)
                    distance_matrix[j, i] = distance_matrix[i, j]
                else:
                    distance_matrix[i, j] = 1.0
                    distance_matrix[j, i] = 1.0
    
    return sample_ids, distance_matrix

def main():
    print("=== Generating Distance Matrix ===")
    
    # Ensure output directory exists
    os.makedirs('output', exist_ok=True)
    
    # Read PHYLIP file
    sequences, n_seqs, seq_len = read_phylip_file('data/faba_fingerprint.phy')
    
    if sequences is None:
        print("Failed to read PHYLIP file")
        return
    
    # Calculate distance matrix
    sample_ids, distance_matrix = calculate_simple_distance_matrix(sequences)
    
    # Save distance matrix as CSV
    dist_df = pd.DataFrame(distance_matrix, index=sample_ids, columns=sample_ids)
    dist_df.to_csv('output/hamming_distance_matrix.csv')
    print("✓ Hamming distance matrix saved: output/hamming_distance_matrix.csv")
    
    # Create PHYLIP format distance matrix
    with open('output/distance_matrix.phy', 'w') as f:
        f.write(f' {len(sample_ids)}\n')
        for i, sample_id in enumerate(sample_ids):
            formatted_id = sample_id.ljust(10)
            distances = ' '.join([f'{dist:.6f}' for dist in distance_matrix[i]])
            f.write(f'{formatted_id} {distances}\n')
    
    print("✓ PHYLIP format distance matrix saved: output/distance_matrix.phy")
    
    # Print summary statistics
    print(f"\nDistance Matrix Summary:")
    print(f"  Minimum distance: {distance_matrix.min():.4f}")
    print(f"  Maximum distance: {distance_matrix.max():.4f}")
    print(f"  Mean distance: {distance_matrix.mean():.4f}")
    print(f"  Standard deviation: {distance_matrix.std():.4f}")

if __name__ == "__main__":
    main()