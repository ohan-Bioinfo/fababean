# Faba Bean SNP QC Pipeline Summary

## Pipeline Parameters
- **Input VCF:** pop.vcf.gz (original)
- **Chromosomes:** chr1, chr2, chr3, chr4, chr5, chr6 (standardized names)
- **Samples:** 21 (constant throughout)
- **Threads:** 55

## Filtering Steps and Results

| Stage | Parameters | Samples | SNPs | Description |
|-------|------------|---------|------|-------------|
| 01_Raw | VCF→PLINK | 21 | 550755 | Original VCF, chr1-chr6 only, SNPs only |
| 02_GENO | --geno 0.05 | 21 | 33077 | Remove variants with >5% missingness |
| 02_MAF | --maf 0.05 | 21 | 33077 | Remove variants with MAF <5% |
| 03_LD_Prune | --indep-pairwise 50 5 0.5 | 21 | 17170 | LD pruning |

## Retention Rates
- **Overall retention:** 3% of raw SNPs
- **GENO filter removed:** 93% of variants
- **LD pruning removed:** 48% of variants

## Final Dataset
- **Total high-quality SNPs:** 17170
- **Samples:** 21
- **Chromosomes:** 6 (chr1-chr6)
