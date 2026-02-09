# scripts/select_top_snps.py
import pandas as pd

def select_top_snps():
    # Read PIC summary
    df_pic = pd.read_csv('output/pic_summary_all.csv')
    
    # Select top 150 SNPs
    top_150 = df_pic.head(150).copy()
    
    # Rename SNPs as SNP1...SNP150
    top_150['New_SNP_ID'] = [f'SNP{i+1:03d}' for i in range(len(top_150))]
    
    # Save the selection
    top_150[['SNP', 'New_SNP_ID', 'CHR', 'A1', 'A2', 'MAF', 'PIC']].to_csv(
        'output/top_150_snps_selection.csv', index=False
    )
    
    # Save SNP list for PLINK
    top_150['SNP'].to_csv('data/top_150_snps_list.txt', 
                         index=False, header=False)
    
    print(f"Selected top 150 SNPs with PIC range: {top_150['PIC'].min():.3f} - {top_150['PIC'].max():.3f}")
    return top_150

if __name__ == "__main__":
    select_top_snps()
