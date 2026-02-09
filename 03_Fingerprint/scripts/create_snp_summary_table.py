import pandas as pd
import numpy as np
import os
from collections import Counter

def fix_pic_calculation():
    """Fix the PIC calculation in the existing diversity file"""
    print("=== FIXING PIC CALCULATION ===")
    
    # Read the existing diversity file
    div_file = 'output/faba_genetic_diversity.csv'
    div_data = pd.read_csv(div_file)
    
    print(f"Original PIC range: {div_data['PIC'].min():.4f} - {div_data['PIC'].max():.4f}")
    
    # Recalculate PIC correctly
    # PIC = 1 - (p² + q²) - 2p²q²
    # For MAF = m, p = m, q = 1-m (if m is the minor allele frequency)
    div_data['PIC_corrected'] = 1 - (div_data['MAF']**2 + (1-div_data['MAF'])**2) - 2 * (div_data['MAF']**2 * (1-div_data['MAF'])**2)
    
    print(f"Corrected PIC range: {div_data['PIC_corrected'].min():.4f} - {div_data['PIC_corrected'].max():.4f}")
    
    # Replace the PIC column
    div_data['PIC'] = div_data['PIC_corrected']
    div_data = div_data.drop('PIC_corrected', axis=1)
    
    # Save corrected diversity file
    div_data.to_csv('output/faba_genetic_diversity_corrected.csv', index=False)
    print("✓ Corrected diversity file saved")
    
    return div_data

def create_final_snp_table():
    """Create the final SNP table using corrected data"""
    print("\n=== CREATING FINAL SNP TABLE ===")
    
    # Read BIM data
    bim_file = 'data/Faba_high_quality.bim'
    bim_data = pd.read_csv(bim_file, 
                          sep='\t', 
                          header=None,
                          names=['CHR', 'SNP', 'cM', 'POS', 'A1', 'A2'])
    
    # Read top SNPs list
    with open('data/top_150_snps_list.txt', 'r') as f:
        top_snps = [line.strip() for line in f if line.strip()]
    
    # Filter for top SNPs
    top_snps_data = bim_data[bim_data['SNP'].isin(top_snps)].copy()
    print(f"✓ Found {len(top_snps_data)} top SNPs in BIM file")
    
    # Use the corrected diversity data
    corrected_div_file = 'output/faba_genetic_diversity_corrected.csv'
    if os.path.exists(corrected_div_file):
        div_data = pd.read_csv(corrected_div_file)
        print(f"✓ Using corrected diversity data for {len(div_data)} SNPs")
        
        # Check SNP name matching
        print(f"BIM SNP format example: {top_snps_data['SNP'].iloc[0]}")
        print(f"Diversity SNP format example: {div_data['Marker'].iloc[0]}")
        
        # Merge using the Marker column from diversity file
        merged_data = top_snps_data.merge(
            div_data, 
            left_on='SNP', 
            right_on='Marker', 
            how='left'
        )
        
        print(f"✓ Merged {len(merged_data)} SNPs")
        
        # Check for any missing values after merge
        missing_pic = merged_data['PIC'].isna().sum()
        missing_gd = merged_data['GeneDiversity'].isna().sum()
        missing_het = merged_data['Heterozygosity'].isna().sum()
        
        print(f"Missing after merge - PIC: {missing_pic}, GeneDiversity: {missing_gd}, Heterozygosity: {missing_het}")
        
    else:
        print("❌ Corrected diversity file not found")
        return
    
    # Create formatted SNP name
    merged_data['SNP_Formatted'] = merged_data.apply(
        lambda row: f"chr{row['CHR']}:{row['POS']}:{row['A1']}:{row['A2']}", axis=1
    )
    
    # Assign new SNP IDs
    merged_data = merged_data.reset_index(drop=True)
    merged_data['New_SNP_ID'] = [f'SNP{str(i+1).zfill(3)}' for i in range(len(merged_data))]
    
    # Create final table with correct columns
    final_columns = [
        'SNP_Formatted', 'New_SNP_ID', 'CHR', 'A1', 'A2', 
        'MAF', 'PIC', 'GeneDiversity', 'Heterozygosity'
    ]
    
    snp_table = merged_data[final_columns].copy()
    
    # Rename columns to match requested format
    snp_table = snp_table.rename(columns={
        'SNP_Formatted': 'SNP',
        'GeneDiversity': 'Genetic_diversity'
    })
    
    # Round to 4 decimal places
    numeric_cols = ['MAF', 'PIC', 'Genetic_diversity', 'Heterozygosity']
    for col in numeric_cols:
        snp_table[col] = snp_table[col].round(4)
    
    # Save final table
    output_file = 'output/SNP_summary_table_final.csv'
    snp_table.to_csv(output_file, index=False)
    
    print(f"\n✅ FINAL SNP TABLE CREATED!")
    print(f"📊 Table contains {len(snp_table)} SNPs")
    print(f"💾 Saved as: {output_file}")
    
    # Display final statistics
    print(f"\n📈 FINAL STATISTICS:")
    for col in numeric_cols:
        if col in snp_table.columns:
            values = snp_table[col]
            print(f"  {col}: {values.min():.4f} - {values.max():.4f} (mean: {values.mean():.4f})")
    
    print(f"\n📋 SAMPLE OF FINAL TABLE:")
    print(snp_table.head(10).to_string(index=False))
    
    return snp_table

def verify_snp_name_matching():
    """Verify that SNP names match between files"""
    print("\n=== VERIFYING SNP NAME MATCHING ===")
    
    # Read all files and compare SNP names
    bim_file = 'data/Faba_high_quality.bim'
    bim_data = pd.read_csv(bim_file, 
                          sep='\t', 
                          header=None,
                          names=['CHR', 'SNP', 'cM', 'POS', 'A1', 'A2'])
    
    # Read top SNPs list
    with open('data/top_150_snps_list.txt', 'r') as f:
        top_snps = [line.strip() for line in f if line.strip()]
    
    # Read diversity file
    div_file = 'output/faba_genetic_diversity.csv'
    div_data = pd.read_csv(div_file)
    
    print(f"BIM file SNPs: {len(bim_data)}")
    print(f"Top SNPs list: {len(top_snps)}")
    print(f"Diversity file SNPs: {len(div_data)}")
    
    # Check overlap
    bim_top = bim_data[bim_data['SNP'].isin(top_snps)]
    div_top = div_data[div_data['Marker'].isin(top_snps)]
    
    print(f"Top SNPs in BIM: {len(bim_top)}")
    print(f"Top SNPs in Diversity: {len(div_top)}")
    
    # Check a few examples
    print(f"\n🔍 SAMPLE SNP NAMES:")
    print(f"BIM: {bim_top['SNP'].iloc[0] if len(bim_top) > 0 else 'N/A'}")
    print(f"Diversity: {div_top['Marker'].iloc[0] if len(div_top) > 0 else 'N/A'}")
    
    # Check if they match exactly
    if len(bim_top) > 0 and len(div_top) > 0:
        same = bim_top['SNP'].iloc[0] == div_top['Marker'].iloc[0]
        print(f"First SNP matches: {same}")
        
        if not same:
            print("❌ SNP names don't match! Need to fix naming convention.")

def create_table_direct_from_diversity():
    """Create table directly from diversity file if merge issues persist"""
    print("\n=== CREATING TABLE DIRECTLY FROM DIVERSITY FILE ===")
    
    # Read the corrected diversity file
    div_file = 'output/faba_genetic_diversity_corrected.csv'
    div_data = pd.read_csv(div_file)
    
    # Extract chromosome and position from Marker
    div_data[['chr_part', 'pos_part']] = div_data['Marker'].str.split(':', expand=True)[[0, 1]]
    div_data['CHR'] = div_data['chr_part'].str.replace('chr', '')
    div_data['POS'] = div_data['pos_part'].astype(int)
    
    # Extract alleles (this is trickier - we need to parse the alleles column)
    div_data[['A1', 'A2']] = div_data['alleles'].str.split('/', expand=True)
    
    # Assign new SNP IDs
    div_data = div_data.reset_index(drop=True)
    div_data['New_SNP_ID'] = [f'SNP{str(i+1).zfill(3)}' for i in range(len(div_data))]
    
    # Create final table
    final_columns = [
        'Marker', 'New_SNP_ID', 'CHR', 'A1', 'A2', 
        'MAF', 'PIC', 'GeneDiversity', 'Heterozygosity'
    ]
    
    snp_table = div_data[final_columns].copy()
    
    # Rename columns
    snp_table = snp_table.rename(columns={
        'Marker': 'SNP',
        'GeneDiversity': 'Genetic_diversity'
    })
    
    # Round to 4 decimal places
    numeric_cols = ['MAF', 'PIC', 'Genetic_diversity', 'Heterozygosity']
    for col in numeric_cols:
        snp_table[col] = snp_table[col].round(4)
    
    # Save this version
    output_file = 'output/SNP_summary_table_direct.csv'
    snp_table.to_csv(output_file, index=False)
    
    print(f"✅ DIRECT SNP TABLE CREATED!")
    print(f"📊 Table contains {len(snp_table)} SNPs")
    print(f"💾 Saved as: {output_file}")
    
    # Display statistics
    print(f"\n📈 DIRECT TABLE STATISTICS:")
    for col in numeric_cols:
        if col in snp_table.columns:
            values = snp_table[col]
            print(f"  {col}: {values.min():.4f} - {values.max():.4f} (mean: {values.mean():.4f})")
    
    print(f"\n📋 SAMPLE OF DIRECT TABLE:")
    print(snp_table.head(10).to_string(index=False))
    
    return snp_table

def fix_plink_hwe_issue():
    """Fix the PLINK HWE column name issue"""
    print("\n=== FIXING PLINK HWE ISSUE ===")
    
    hwe_file = 'output/faba_150_hwe.hwe'
    if os.path.exists(hwe_file):
        hwe_data = pd.read_csv(hwe_file, sep='\s+')
        print(f"Columns in HWE file: {hwe_data.columns.tolist()}")
        
        # PLINK HWE file typically has columns: CHR, SNP, TEST, A1, A2, GENO, O(HET), E(HET), P
        # We need to calculate allele frequencies from the genotype counts
        if 'O(HET)' in hwe_data.columns and 'GENO' in hwe_data.columns:
            # Parse GENO column which has format like: 4/7/10 (hom_ref/het/hom_alt)
            hwe_data[['hom_ref', 'het', 'hom_alt']] = hwe_data['GENO'].str.split('/', expand=True).astype(int)
            
            # Calculate total individuals and alleles
            hwe_data['N'] = hwe_data['hom_ref'] + hwe_data['het'] + hwe_data['hom_alt']
            hwe_data['ref_alleles'] = hwe_data['hom_ref'] * 2 + hwe_data['het']
            hwe_data['alt_alleles'] = hwe_data['hom_alt'] * 2 + hwe_data['het']
            hwe_data['total_alleles'] = hwe_data['ref_alleles'] + hwe_data['alt_alleles']
            
            # Calculate allele frequencies
            hwe_data['P'] = hwe_data['ref_alleles'] / hwe_data['total_alleles']
            hwe_data['Q'] = hwe_data['alt_alleles'] / hwe_data['total_alleles']
            
            # Calculate PIC and Gene Diversity
            hwe_data['PIC'] = 1 - (hwe_data['P']**2 + hwe_data['Q']**2) - 2 * (hwe_data['P']**2 * hwe_data['Q']**2)
            hwe_data['GeneDiversity'] = 2 * hwe_data['P'] * hwe_data['Q']
            
            print(f"PIC from HWE: {hwe_data['PIC'].min():.4f} - {hwe_data['PIC'].max():.4f}")
            print(f"GeneDiversity from HWE: {hwe_data['GeneDiversity'].min():.4f} - {hwe_data['GeneDiversity'].max():.4f}")
            
            # Save corrected HWE data
            hwe_data.to_csv('output/faba_150_hwe_corrected.csv', index=False)
            print("✓ Corrected HWE data saved")

if __name__ == "__main__":
    # First, verify SNP name matching
    verify_snp_name_matching()
    
    # Fix the PIC calculation in the diversity file
    corrected_div = fix_pic_calculation()
    
    # Try to create the final table with proper merging
    final_table = create_final_snp_table()
    
    # If that still has issues, create table directly from diversity file
    if final_table is None or final_table['PIC'].isna().all():
        print("\n⚠️  Merge issues detected, creating table directly from diversity file...")
        direct_table = create_table_direct_from_diversity()
    
    # Fix the PLINK HWE issue
    fix_plink_hwe_issue()