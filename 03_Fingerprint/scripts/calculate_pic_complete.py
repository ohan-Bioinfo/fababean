# scripts/calculate_pic_complete.py
import pandas as pd
import numpy as np
import subprocess
import os

def calculate_pic(p, q):
    """Calculate Polymorphic Information Content for biallelic SNP"""
    return 1 - (p**2 + q**2) - 2*(p**2)*(q**2)

def main():
    # Read allele frequencies
    frq_file = "data/Faba_high_quality.frq"
    
    if not os.path.exists(frq_file):
        print("Generating allele frequencies with PLINK...")
        subprocess.run([
            "plink", "--bfile", "data/Faba_high_quality", 
            "--freq", "--out", "data/Faba_high_quality"
        ])
    
    # Read frequency file
    df_frq = pd.read_csv(frq_file, delim_whitespace=True)
    
    # Calculate PIC
    df_frq['PIC'] = df_frq.apply(
        lambda row: calculate_pic(row['MAF'], 1 - row['MAF']), 
        axis=1
    )
    
    # Sort by PIC (descending)
    df_frq_sorted = df_frq.sort_values('PIC', ascending=False)
    
    # Save PIC summary
    df_frq_sorted[['CHR', 'SNP', 'A1', 'A2', 'MAF', 'PIC']].to_csv(
        'output/pic_summary_all.csv', index=False
    )
    
    print(f"Total SNPs with PIC calculated: {len(df_frq_sorted)}")
    print(f"Top 5 SNPs by PIC:")
    print(df_frq_sorted.head()[['SNP', 'MAF', 'PIC']])
    
    return df_frq_sorted

if __name__ == "__main__":
    df_pic = main()
