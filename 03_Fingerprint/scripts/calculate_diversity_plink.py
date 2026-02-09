import pandas as pd
import numpy as np
import subprocess
import os

def calculate_diversity_from_plink():
    """Calculate diversity statistics using PLINK commands"""
    
    # Use PLINK to calculate allele frequencies
    plink_cmd = [
        "plink",
        "--bfile", "../data/Faba_high_quality",
        "--freq",
        "--out", "../output/faba_allele_freq"
    ]
    
    subprocess.run(plink_cmd)
    
    # Use PLINK to calculate heterozygosity
    plink_cmd = [
        "plink", 
        "--bfile", "../data/Faba_high_quality",
        "--het",
        "--out", "../output/faba_heterozygosity"
    ]
    
    subprocess.run(plink_cmd)
    
    # Read and process results
    freq_df = pd.read_csv("../output/faba_allele_freq.frq", delim_whitespace=True)
    het_df = pd.read_csv("../output/faba_heterozygosity.het", delim_whitespace=True)
    
    # Calculate diversity statistics
    results = []
    for _, row in freq_df.iterrows():
        maf = min(row['MAF'], 1 - row['MAF'])
        gene_diversity = 2 * maf * (1 - maf)
        pic = 1 - (maf**2 + (1-maf)**2) - 2 * (maf**2 * (1-maf)**2)
        
        results.append({
            'Marker': row['SNP'],
            'chr': row['CHR'],
            'position': row['POS'],
            'alleles': f"{row['A1']}/{row['A2']}",
            'PIC': round(pic, 4),
            'MAF': round(maf, 4),
            'GeneDiversity': round(gene_diversity, 4),
            'Heterozygosity': 'NA'  # Would need individual level data
        })
    
    df_results = pd.DataFrame(results)
    df_results.to_csv("../output/faba_genetic_diversity_plink.csv", index=False)
    
    print(f"Genetic diversity analysis of {len(df_results)} markers")
    print("="*80)
    print(df_results.to_string(index=False))
    
    return df_results

if __name__ == "__main__":
    calculate_diversity_from_plink()
