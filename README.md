# Faba Bean GBS Analysis (iter-3)

This repository contains a stage-organized Faba bean (Vicia faba) GBS analysis workflow and results, prepared for GitHub publication and reproducible reruns.

## Project Layout

- `00_GBS/`: input VCF summaries and variant overview plots.
- `01_Raw/`: PLINK conversion, QC filtering, LD pruning, and QC plots/tables.
- `02_Diversity/`: PCA, IBS/IBD, kinship, and ADMIXTURE analysis.
- `03_Fingerprint/`: top-SNP fingerprint panel generation and genotype heatmaps.
- `04_PhylogeneticTree/`: NJ/ML phylogenetic analyses and tree visualizations.
- `workflow/`: reproducible stage runner scripts.
- `envs/environment.yml`: conda environment for tools and libraries.
- `docs/STAGE_INVENTORY.md`: stage-by-stage inputs/outputs summary.

## Reproducible Execution

Run from repository root.

```bash
conda env create -f envs/environment.yml
conda activate fababean-gbs
```

Run all stages:

```bash
bash workflow/run_all.sh
```

Run one stage:

```bash
bash workflow/stage_01_qc.sh
bash workflow/stage_02_diversity.sh
bash workflow/stage_03_fingerprint.sh
bash workflow/stage_04_phylogeny.sh
```

Optional settings:

```bash
THREADS=16 bash workflow/run_all.sh
DRY_RUN=1 bash workflow/run_all.sh
K_MIN=2 K_MAX=10 ADMIXTURE_SEED=43 bash workflow/stage_02_diversity.sh
```

## Notes

- `workflow/stage_02_diversity.sh` regenerates `02_Diversity/Admixture/cv_summary.txt` from `log_K*.out`.
- `02_Diversity/Admixture/plot_cv_and_admixture.R` appears to be a shell transcript; use `02_Diversity/Admixture/visualize_cv.R` for plotting.
- `04_PhylogeneticTree/scripts/build_ml_tree.py` is empty; reproducible ML analysis uses the comprehensive runner instead.

## Publish To GitHub

```bash
git init
git add .
git commit -m "Add reproducible Faba bean GBS analysis workflow"
git branch -M main
git remote add origin <your-github-repo-url>
git push -u origin main
```

If repository size grows, enable Git LFS for large binary outputs (VCF/PLINK/figures).
