import pandas as pd
import numpy as np
from collections import Counter
import subprocess
import os

def calculate_heterozygosity(genotypes):
    """Calculate observed heterozygosity"""
    total = len(genotypes)
    if total == 0:
        return 0
    het_count = sum(1 for gt in genotypes if len(set(gt)) > 1)
    return het_count / total

def calculate_expected_heterozygosity(allele_freq):
    """Calculate expected heterozygosity (Gene Diversity)"""
    if sum(allele_freq.values()) == 0:
        return 0
    return 1 - sum((freq/sum(allele_freq.values()))**2 for freq in allele_freq.values())

def calculate_pic(allele_freq):
    """Calculate Polymorphism Information Content"""
    freqs = [f/sum(allele_freq.values()) for f in allele_freq.values() if f > 0]
    if len(freqs) < 2:
        return 0
    sum_sq = sum(f**2 for f in freqs)
    return 1 - sum_sq - sum(2 * f1 * f2 for i, f1 in enumerate(freqs) for f2 in freqs[i+1:])

def calculate_maf(allele_freq):
    """Calculate Minor Allele Frequency"""
    if not allele_freq:
        return 0
    total = sum(allele_freq.values())
    freqs = [f/total for f in allele_freq.values()]
    return min(freqs)

def get_allele_frequencies(genotypes):
    """Calculate allele frequencies from genotypes"""
    alleles = []
    for gt in genotypes:
        alleles.extend(gt)
    return Counter(alleles)

def analyze_marker_diversity(vcf_file, output_file):
    """Analyze genetic diversity for each marker in VCF file"""
    
    # Read VCF file
    vcf_data = []
    with open(vcf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            vcf_data.append(line.strip().split('\t'))
    
    # Extract sample names (assuming they start from column 9)
    samples = []
    if vcf_data:
        with open(vcf_file, 'r') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    samples = line.strip().split('\t')[9:]
                    break
    
    results = []
    
    for row in vcf_data:
        chrom = row[0]
        pos = row[1]
        marker_id = row[2] if row[2] != '.' else f"{chrom}_{pos}"
        ref_allele = row[3]
        alt_alleles = row[4].split(',')
        
        # All possible alleles
        all_alleles = [ref_allele] + alt_alleles
        
        # Parse genotypes
        genotypes = []
        for sample_gt in row[9:]:
            gt_field = sample_gt.split(':')[0]  # Take GT field (first part)
            if '/' in gt_field:
                alleles_idx = gt_field.split('/')
            elif '|' in gt_field:
                alleles_idx = gt_field.split('|')
            else:
                continue
                
            # Convert to actual alleles
            try:
                genotype = []
                for idx in alleles_idx:
                    if idx != '.':
                        genotype.append(all_alleles[int(idx)])
                if genotype:
                    genotypes.append(genotype)
            except (ValueError, IndexError):
                continue
        
        if not genotypes:
            continue
            
        # Calculate statistics
        allele_freq = get_allele_frequencies(genotypes)
        alleles_str = '/'.join(allele_freq.keys())
        maf = calculate_maf(allele_freq)
        gene_diversity = calculate_expected_heterozygosity(allele_freq)
        heterozygosity = calculate_heterozygosity(genotypes)
        pic = calculate_pic(allele_freq)
        
        results.append({
            'Marker': marker_id,
            'chr': chrom,
            'position': pos,
            'alleles': alleles_str,
            'PIC': round(pic, 4),
            'MAF': round(maf, 4),
            'GeneDiversity': round(gene_diversity, 4),
            'Heterozygosity': round(heterozygosity, 4)
        })
    
    # Create results DataFrame
    df_results = pd.DataFrame(results)
    
    # Save to file
    df_results.to_csv(output_file, index=False)
    print(f"Results saved to {output_file}")
    
    # Print summary
    print(f"\nGenetic diversity analysis of {len(df_results)} markers in {len(samples)} faba bean materials")
    print("="*80)
    print(df_results.to_string(index=False))
    
    return df_results

if __name__ == "__main__":
    vcf_file = "../output/faba_150.vcf"
    output_file = "../output/faba_genetic_diversity.csv"
    
    if os.path.exists(vcf_file):
        analyze_marker_diversity(vcf_file, output_file)
    else:
        print(f"VCF file {vcf_file} not found!")
        print("Available files in output directory:")
        for file in os.listdir("../output"):
            if file.endswith('.vcf') or file.endswith('.csv'):
                print(f"  - {file}")
