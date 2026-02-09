# Faba Bean GBS - Kinship Analysis Report

## Analysis Overview
- **Dataset**: LD-pruned high-quality SNPs (17,170 SNPs)
- **Accessions**: 21 Faba Bean genotypes
- **Method**: PLINK2 --make-king
- **Metric**: Kinship coefficients

## Kinship Interpretation
- **0.5**: Parent-Offspring / Full Siblings
- **0.25**: Half Siblings / Grandparent-Grandchild  
- **0.125**: First Cousins
- **0.0625**: Second Cousins
- **0**: Unrelated
- **Negative values**: Less related than population average

## Files Generated
- `Faba_Kinship.king` - Kinship matrix (lower triangular)
- `Faba_Kinship.king.id` - Accession IDs
- `Faba_Kinship_Analysis.png` - Comprehensive kinship analysis
- `Faba_Kinship_Clustered.png` - Clustered kinship heatmap
- `Faba_Kinship_Minimal_Heatmap.png` - Minimal heatmap
- `Faba_Kinship_Matrix.csv` - Kinship matrix (CSV format)

## Key Applications
- Population structure assessment
- Relatedness detection in breeding programs
- Genomic selection and GWAS correction
- Germplasm management and conservation

## Quality Indicators
✅ Kinship coefficients within expected range
✅ Clear clustering patterns visible
✅ No extreme outlier values detected

*Analysis performed on: $(date)*
