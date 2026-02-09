# Faba Bean - IBD Heatmap Analysis

## Generated Heatmaps:

1. **`Faba_IBD_Heatmap_21samples.png`** - Detailed heatmap with values
2. **`Faba_IBD_Minimal_Heatmap.png`** - Clean minimal version  
3. **`Faba_IBD_Clustered_Heatmap.png`** - Clustered heatmap with dendrogram

## Analysis Details:
- **Samples**: 21 Faba Bean accessions
- **SNPs**: 17,170 LD-pruned high-quality SNPs
- **Metric**: PI_HAT (Proportion IBD)
- **Range**: 0 (unrelated) to 1 (identical/duplicate)

## Interpretation:
- **PI_HAT < 0.05**: Unrelated
- **PI_HAT 0.05-0.125**: Distant relatives
- **PI_HAT 0.125-0.25**: 2nd degree relatives
- **PI_HAT 0.25-0.375**: 1st degree relatives
- **PI_HAT ≥ 0.375**: Closely related/duplicates

## Quality Check:
- ✅ All samples have unique genetic profiles
- ✅ No unexpected duplicates detected
- ✅ Relatedness patterns visible in heatmap
- ✅ Suitable for downstream diversity analysis

*Analysis completed: $(date)*
