# Faba Bean GBS - SNP Density Analysis (5Mb Windows)

## Analysis Overview
- **Dataset**: LD-pruned high-quality SNPs (17,170 total)
- **Window size**: 5 Mb
- **Chromosomes**: chr1L, chr1S, 2, 3, 4, 5, 6

## Files Generated
- `Faba_Bean_SNP_Density_Heatmap_5Mb.png` - Visual heatmap of SNP density
- `snp_density_5mb_windows.csv` - Raw density data for further analysis

## Expected Patterns to Observe

1. **High-Density Regions**: Likely gene-rich regions or recombination hotspots
2. **Low-Density Regions**: Centromeres, telomeres, or repetitive regions
3. **Chromosome Arms**: Often show different density patterns
4. **Distribution**: Should be relatively even across chromosomes

## Interpretation Guide

| Density Range | Interpretation |
|---------------|----------------|
| 0 SNPs/window | Possible assembly gaps or low-complexity regions |
| 1-10 SNPs/window | Low density, possibly repetitive regions |
| 11-50 SNPs/window | Moderate density |
| 51+ SNPs/window | High density, potential gene-rich regions |

## Quality Indicators
- ✅ Even distribution across chromosomes
- ✅ No large contiguous zero-SNP regions (except centromeres)
- ✅ Reasonable density variation (10-100 SNPs/5Mb typical for GBS)

*Note: GBS data typically shows uneven distribution due to restriction enzyme cut sites*
