# scripts/convert_raw_to_phylip.py
import pandas as pd
import numpy as np
import subprocess
import os

def generate_raw_file():
    """Generate PLINK RAW file if it doesn't exist"""
    print("Generating PLINK RAW file...")
    
    # Check if we already have the files
    if os.path.exists('data/faba_fingerprint_fixed.raw'):
        print("✓ RAW file already exists")
        return True
    
    # Run PLINK to generate RAW file
    cmd = [
        'plink',
        '--bfile', 'data/faba_fingerprint',
        '--allow-extra-chr',
        '--recode', 'A',
        '--out', 'data/faba_fingerprint_fixed'
    ]
    
    print(f"Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode == 0:
        print("✓ RAW file generated successfully")
        return True
    else:
        print(f"✗ Failed to generate RAW file: {result.stderr}")
        return False

def convert_to_phylip():
    """Convert RAW file to PHYLIP format"""
    print("Converting RAW file to PHYLIP format...")
    
    # Read RAW file
    try:
        raw = pd.read_csv('data/faba_fingerprint_fixed.raw', sep='\s+')
        print(f"✓ RAW file loaded: {len(raw)} samples")
    except Exception as e:
        print(f"✗ Error reading RAW file: {e}")
        return False
    
    # Extract sample IDs and SNP columns
    sample_ids = raw['IID'].tolist()
    snp_cols = [col for col in raw.columns if col.startswith('SNP')]
    
    print(f"✓ Found {len(snp_cols)} SNPs")
    
    # Create PHYLIP format
    with open('data/faba_fingerprint.phy', 'w') as f:
        # Header: number of samples and number of sites
        f.write(f' {len(sample_ids)} {len(snp_cols)}\n')
        
        # Write each sample's genotypes
        for i, sample_id in enumerate(sample_ids):
            # Format sample ID (max 10 characters for PHYLIP)
            formatted_id = sample_id[:10].ljust(10)
            
            # Convert genotypes to characters (0,1,2) and handle missing data
            genotype_chars = []
            for snp in snp_cols:
                val = raw[snp].iloc[i]
                if pd.isna(val):
                    genotype_chars.append('?')
                else:
                    genotype_chars.append(str(int(val)))
            
            genotype_str = ''.join(genotype_chars)
            f.write(f'{formatted_id} {genotype_str}\n')
    
    print(f"✓ PHYLIP file created: data/faba_fingerprint.phy")
    
    # Verify the file
    with open('data/faba_fingerprint.phy', 'r') as f:
        lines = f.readlines()
        print(f"✓ PHYLIP file verified: {len(lines)} lines")
        print("First 3 lines:")
        for i in range(min(3, len(lines))):
            print(f"  {lines[i].strip()}")
    
    return True

def main():
    print("=== Alternative PHYLIP Conversion ===")
    
    # Ensure data directory exists
    os.makedirs('data', exist_ok=True)
    
    # Generate RAW file
    if generate_raw_file():
        # Convert to PHYLIP
        if convert_to_phylip():
            print("\n✓ Successfully created PHYLIP file via alternative method")
        else:
            print("\n✗ Alternative conversion failed")
    else:
        print("\n✗ Cannot generate RAW file")

if __name__ == "__main__":
    main()