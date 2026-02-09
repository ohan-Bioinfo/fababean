import pandas as pd
import numpy as np
import os
import warnings
warnings.filterwarnings('ignore')

class SNPDataLoader:
    def __init__(self, base_path):
        self.base_path = base_path
        self.data = {}
    
    def load_all_data(self):
        """Load all SNP data files with proper error handling"""
        print("Loading SNP data...")
        
        try:
            # Load missingness data with proper separator
            self.data['sample_missingness'] = pd.read_csv(
                os.path.join(self.base_path, '01_SampleMissingness_QC/Faba_chrOnly_raw.imiss'), 
                sep='\s+', comment='#'
            )
            self.data['snp_missingness'] = pd.read_csv(
                os.path.join(self.base_path, '01_SampleMissingness_QC/Faba_chrOnly_raw.lmiss'), 
                sep='\s+', comment='#'
            )
            print("✓ Missingness data loaded")
            
        except Exception as e:
            print(f"✗ Error loading missingness data: {e}")
            # Try alternative loading method
            try:
                self.data['sample_missingness'] = pd.read_csv(
                    os.path.join(self.base_path, '01_SampleMissingness_QC/Faba_chrOnly_raw.imiss'), 
                    delim_whitespace=True
                )
                self.data['snp_missingness'] = pd.read_csv(
                    os.path.join(self.base_path, '01_SampleMissingness_QC/Faba_chrOnly_raw.lmiss'), 
                    delim_whitespace=True
                )
                print("✓ Missingness data loaded with alternative method")
            except Exception as e2:
                print(f"✗ Alternative loading also failed: {e2}")
        
        try:
            # Load heterozygosity data - check the actual column names
            het_file = os.path.join(self.base_path, '04_Het_QC/Faba_chrOnly_het.het')
            if os.path.exists(het_file):
                het_df = pd.read_csv(het_file, sep='\s+', comment='#')
                print(f"Heterozygosity file columns: {het_df.columns.tolist()}")
                self.data['heterozygosity'] = het_df
                print("✓ Heterozygosity data loaded")
        except Exception as e:
            print(f"✗ Error loading heterozygosity data: {e}")
        
        try:
            # Load BIM files for SNP counts
            self.data['raw_bim'] = pd.read_csv(
                os.path.join(self.base_path, 'Faba_chrOnly_raw.bim'), 
                sep='\s+', header=None, 
                names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2']
            )
            self.data['geno05_bim'] = pd.read_csv(
                os.path.join(self.base_path, '02_SNP_Filter/Faba_chrOnly_geno05.bim'), 
                sep='\s+', header=None, 
                names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2']
            )
            self.data['maf05_bim'] = pd.read_csv(
                os.path.join(self.base_path, '02_SNP_Filter/Faba_chrOnly_geno05_maf05.bim'), 
                sep='\s+', header=None, 
                names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2']
            )
            self.data['pruned_bim'] = pd.read_csv(
                os.path.join(self.base_path, '03_LD_Prune/Faba_chrOnly_pruned.bim'), 
                sep='\s+', header=None, 
                names=['chr', 'snp_id', 'cm', 'pos', 'a1', 'a2']
            )
            print("✓ BIM files loaded")
            
        except Exception as e:
            print(f"✗ Error loading BIM files: {e}")
        
        # Print summary of loaded data
        print("\n=== Data Loading Summary ===")
        for key, value in self.data.items():
            if hasattr(value, 'shape'):
                print(f"{key}: {value.shape}")
            else:
                print(f"{key}: Loaded")
        
        return self.data

    def get_summary_stats(self):
        """Generate basic summary statistics"""
        if not self.data:
            print("No data loaded. Call load_all_data() first.")
            return None
            
        summary = {
            'Total_Raw_SNPs': len(self.data.get('raw_bim', [])),
            'Total_Samples': len(self.data.get('sample_missingness', [])),
            'SNPs_after_Geno_Filter': len(self.data.get('geno05_bim', [])),
            'SNPs_after_MAF_Filter': len(self.data.get('maf05_bim', [])),
            'SNPs_after_LD_Pruning': len(self.data.get('pruned_bim', [])),
        }
        
        if 'sample_missingness' in self.data:
            summary['Mean_Sample_Missingness'] = self.data['sample_missingness']['F_MISS'].mean()
        
        return summary

# Test the data loader
if __name__ == "__main__":
    loader = SNPDataLoader(".")
    data = loader.load_all_data()
    stats = loader.get_summary_stats()
    print("\nSummary Statistics:")
    for key, value in stats.items():
        print(f"{key}: {value}")
