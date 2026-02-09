# Faba Bean GBS - IBD Analysis Report

## Analysis Overview
- **Dataset**: LD-pruned high-quality SNPs (17,170 SNPs)
- **Samples**: 21
- **Method**: PLINK --genome (IBD estimation)
- **Pairs analyzed**: 210 (n*(n-1)/2 for n=21 samples)

## Relatedness Interpretation
- **PI_HAT < 0.05**: Unrelated
- **PI_HAT 0.05-0.125**: 3rd degree relatives
- **PI_HAT 0.125-0.25**: 2nd degree relatives  
- **PI_HAT 0.25-0.375**: 1st degree relatives
- **PI_HAT ≥ 0.375**: Duplicates/MZ twins

## Quality Control
- ✅ Check for unexpected relatedness
- ✅ Identify potential duplicates
- ✅ Assess overall population structure
- ✅ Validate sample relationships

## Files Generated
- `Faba_IBD.genome` - Raw IBD estimates from PLINK
- `Faba_IBD_Analysis.png` - Comprehensive visualization
- `IBD_summary_stats.txt` - Statistical summary
- `IBD_analysis_report.md` - This report

## Next Steps
1. Review closely related pairs
2. Consider removing duplicates if found
3. Proceed to PCA and kinship analysis
4. Use IBD information for population structure correction

*Analysis performed on: $(date)*
