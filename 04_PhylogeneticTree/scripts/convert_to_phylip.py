# scripts/convert_to_phylip_fixed.py
import pandas as pd
import numpy as np
import subprocess
import os

def fix_plink_files():
    """Fix PLINK files and generate necessary formats"""
    
    print("Step 1: Fixing PLINK files and generating required formats...")
    
    # Generate PED file with proper chromosome handling
    print("Generating PED file...")
    result = subprocess.run([
        'plink', '--bfile', 'data/faba_fingerprint',
        '--allow-extra-chr',
        '--recode',
        '--out', 'data/faba_fingerprint_fixed'
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Error generating PED file:")
        print(result.stderr)
        return False
    
    # Generate RAW file with proper chromosome handling  
    print("Generating RAW file...")
    result = subprocess.run([
        'plink', '--bfile', 'data/faba_fingerprint', 
        '--allow-extra-chr',
        '--recode', 'A',
        '--out', 'data/faba_fingerprint_fixed'
    ], capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Error generating RAW file:")
        print(result.stderr)
        return False
    
    return True

def plink_to_phylip():
    """Convert PLINK data to PHYLIP format for phylogenetic analysis"""
    
    print("Step 2: Converting to PHYLIP format...")
    
    # Read the fixed BIM file
    bim = pd.read_csv('data/faba_fingerprint.bim', 
                     sep='\t', header=None,
                     names=['CHR', 'SNP', 'cM', 'POS', 'A1', 'A2'])
    
    # Read the FAM file to get sample information
    fam = pd.read_csv('data/faba_fingerprint.fam',
                     sep='\t', header=None,
                     names=['FID', 'IID', 'Father', 'Mother', 'Sex', 'Phenotype'])
    
    # Read the fixed PED file
    try:
        ped = pd.read_csv('data/faba_fingerprint_fixed.ped', 
                         sep='\t', header=None)
        print("Successfully read PED file")
    except Exception as e:
        print(f"Error reading PED file: {e}")
        return None
    
    # Read the fixed RAW file for numerical genotypes
    try:
        raw = pd.read_csv('data/faba_fingerprint_fixed.raw', sep='\s+')
        print("Successfully read RAW file")
    except Exception as e:
        print(f"Error reading RAW file: {e}")
        return None
    
    # Extract sample IDs
    sample_ids = fam['IID'].tolist()
    
    # Extract numerical genotypes (columns starting with SNP)
    snp_cols = [col for col in raw.columns if col.startswith('SNP')]
    numeric_genotypes = raw[snp_cols]
    
    print(f"Found {len(sample_ids)} samples and {len(snp_cols)} SNPs")
    
    # Create PHYLIP format
    n_samples = len(raw)
    n_snps = len(snp_cols)
    
    # Create PHYLIP file
    with open('data/faba_fingerprint.phy', 'w') as f:
        # Header: number of samples and number of sites
        f.write(f'{n_samples} {n_snps}\n')
        
        # Write each sample's genotypes
        for i, sample_id in enumerate(raw['IID']):
            # Format sample ID (max 10 characters for PHYLIP)
            formatted_id = sample_id[:10].ljust(10)
            
            # Get genotype string (convert to 0,1,2 but as characters)
            genotype_str = ''.join(numeric_genotypes.iloc[i].astype(int).astype(str))
            
            f.write(f'{formatted_id} {genotype_str}\n')
    
    print(f"PHYLIP file created: {n_samples} samples, {n_snps} SNPs")
    
    # Also create NEXUS format
    create_nexus_format(raw, sample_ids, snp_cols, numeric_genotypes)
    
    return raw, sample_ids, snp_cols

def create_nexus_format(raw, sample_ids, snp_cols, numeric_genotypes):
    """Create NEXUS format file for phylogenetic analysis"""
    
    with open('data/faba_fingerprint.nex', 'w') as f:
        f.write('#NEXUS\n\n')
        f.write('BEGIN DATA;\n')
        f.write(f'    DIMENSIONS NTAX={len(sample_ids)} NCHAR={len(snp_cols)};\n')
        f.write('    FORMAT DATATYPE=STANDARD SYMBOLS="012" MISSING=?;\n')
        f.write('    MATRIX\n')
        
        for i, sample_id in enumerate(sample_ids):
            genotype_str = ''.join(numeric_genotypes.iloc[i].astype(int).astype(str))
            f.write(f'    {sample_id[:20].ljust(20)} {genotype_str}\n')
        
        f.write('    ;\n')
        f.write('END;\n\n')
    
    print("NEXUS file created")

def create_simple_phylip():
    """Create a simple PHYLIP file as fallback"""
    print("Creating simple PHYLIP format as fallback...")
    
    # Read the RAW file directly
    try:
        raw = pd.read_csv('data/faba_fingerprint_fixed.raw', sep='\s+')
    except:
        print("Cannot read RAW file. Please check if PLINK ran successfully.")
        return
    
    sample_ids = raw['IID'].tolist()
    snp_cols = [col for col in raw.columns if col.startswith('SNP')]
    
    # Simple PHYLIP format
    with open('data/simple_fingerprint.phy', 'w') as f:
        f.write(f'{len(sample_ids)} {len(snp_cols)}\n')
        for i, sample_id in enumerate(sample_ids):
            formatted_id = sample_id[:10].ljust(10)
            genotypes = ''.join(raw[snp_cols].iloc[i].astype(int).astype(str))
            f.write(f'{formatted_id} {genotypes}\n')
    
    print("Simple PHYLIP file created as backup")

if __name__ == "__main__":
    print("=== Starting Data Conversion ===")
    
    # First fix the PLINK files
    success = fix_plink_files()
    
    if success:
        # Then convert to phylogenetic formats
        result = plink_to_phylip()
        if result is None:
            print("Trying fallback method...")
            create_simple_phylip()
    else:
        print("Failed to fix PLINK files. Trying alternative approach...")
        create_simple_phylip()
    
    print("=== Data Conversion Complete ===")