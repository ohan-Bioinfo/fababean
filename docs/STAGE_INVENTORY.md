# Stage Inventory

## 00_GBS

- Purpose: source variant summary and initial descriptive plots.
- Key inputs: `00_GBS/pop.vcf.gz`, archived stats in `00_GBS/arc/`.
- Key scripts: `00_GBS/vs.py`.
- Key outputs: `00_GBS/arc/Faba_variant_summary.tsv`, `00_GBS/arc/Faba_variant_summary.png`.

## 01_Raw

- Purpose: convert VCF to PLINK, QC filtering, heterozygosity, LD pruning.
- Key inputs: `00_GBS/pop.vcf.gz`.
- Reproducible runner: `workflow/stage_01_qc.sh`.
- Key outputs:
  - `01_Raw/Faba_chrOnly_raw.{bed,bim,fam}`
  - `01_Raw/02.2__SNP_Filter/Faba_chrOnly_geno05_maf05.{bed,bim,fam}`
  - `01_Raw/03_LD_Prune/Faba_chrOnly_pruned.{bed,bim,fam}`

## 02_Diversity

- Purpose: population structure and relationship analysis (PCA, IBS/IBD, kinship, ADMIXTURE).
- Main input: `01_Raw/03_LD_Prune/Faba_chrOnly_pruned.{bed,bim,fam}`.
- Reproducible runner: `workflow/stage_02_diversity.sh`.
- Key outputs:
  - `02_Diversity/PCA/Faba_PCA.eigenvec`
  - `02_Diversity/IBS/Faba_IBS.mibs`
  - `02_Diversity/IBD/Faba_IBD.genome`
  - `02_Diversity/Kinship/Faba_Kinship.king`
  - `02_Diversity/Admixture/Faba_chrOnly_pruned_numeric.*`
  - `02_Diversity/Admixture/cv_summary.txt`

## 03_Fingerprint

- Purpose: top-SNP fingerprint panel generation and genotype heatmap output.
- Main input: `01_Raw/03_LD_Prune/Faba_chrOnly_pruned.{bed,bim,fam}`.
- Reproducible runner: `workflow/stage_03_fingerprint.sh`.
- Key outputs:
  - `03_Fingerprint/data/faba_fingerprint.{bed,bim,fam}`
  - `03_Fingerprint/output/top_150_snps_selection.csv`
  - `03_Fingerprint/plots/faba_fingerprint_heatmap_categorical.png`

## 04_PhylogeneticTree

- Purpose: phylogenetic reconstruction and tree visualization.
- Main input: `03_Fingerprint/data/faba_fingerprint.*`.
- Reproducible runner: `workflow/stage_04_phylogeny.sh`.
- Key outputs:
  - `04_PhylogeneticTree/data/faba_fingerprint.phy`
  - `04_PhylogeneticTree/output/*.newick`
  - `04_PhylogeneticTree/plots/*.png`
