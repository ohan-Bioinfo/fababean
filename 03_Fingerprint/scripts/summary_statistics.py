import pandas as pd
import numpy as np

def generate_summary_report(diversity_file):
    """Generate summary statistics report"""
    
    df = pd.read_csv(diversity_file)
    
    print("FABA BEAN GENETIC DIVERSITY SUMMARY")
    print("="*50)
    print(f"Total markers analyzed: {len(df)}")
    print(f"Total chromosomes covered: {df['chr'].nunique()}")
    print(f"Average PIC: {df['PIC'].mean():.4f}")
    print(f"Average MAF: {df['MAF'].mean():.4f}")
    print(f"Average Gene Diversity: {df['GeneDiversity'].mean():.4f}")
    print(f"Average Heterozygosity: {df['Heterozygosity'].mean():.4f}")
    print("\n" + "="*50)
    
    # Save summary
    summary = {
        'Statistic': ['Total Markers', 'Chromosomes', 'Mean PIC', 'Mean MAF', 
                     'Mean Gene Diversity', 'Mean Heterozygosity'],
        'Value': [len(df), df['chr'].nunique(), round(df['PIC'].mean(), 4),
                 round(df['MAF'].mean(), 4), round(df['GeneDiversity'].mean(), 4),
                 round(df['Heterozygosity'].mean(), 4)]
    }
    
    summary_df = pd.DataFrame(summary)
    summary_df.to_csv(diversity_file.replace('.csv', '_summary.csv'), index=False)

# Run the scripts in order:
print("1. Running genetic diversity analysis...")
os.system("python calculate_faba_diversity.py")

print("\n2. Generating summary statistics...")
generate_summary_report("../output/faba_genetic_diversity.csv")
