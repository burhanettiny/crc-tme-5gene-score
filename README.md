# A 5-gene tumour microenvironment prognostic score in colorectal cancer

Reproducible R code for:

> Yalçınkaya B. *Development and Validation of a 5-Gene Tumor Microenvironment Prognostic Score in Colorectal Cancer.* Revista da Associação Médica Brasileira (manuscript RAMB-2026-0966, under revision).

The 5-gene score (IL6, CXCL8, PTGS2, TIMP1, MMP3) was derived with LASSO-Cox in TCGA-COAD (n = 436) from tumour-microenvironment hub genes of a discovery PPI network. It was validated in GSE39582 (n = 562 primary tumours).

```
risk score = 0.3558·IL6 − 0.0656·CXCL8 − 0.1132·PTGS2 + 0.2495·TIMP1 − 0.1623·MMP3
```

The coefficients are in `data/processed/tcga_lasso_coefficients.csv`.

## Repository structure

| Path | Content |
|---|---|
| `R/00_helpers.R` | Shared data loading and helper functions (sourced by scripts 08–13) |
| `R/01_cohort1_DEG_PPI.R` | Discovery Cohort-1 (GSE41328, GSE37364, GSE15960): normalization, ComBat, limma, GO/KEGG, STRING PPI |
| `R/02_cohort2_DEG_PPI.R` | Discovery Cohort-2 (GSE21510, GSE32323, GSE110224), including detection of the 17 duplicate patients |
| `R/03_tcga_cox_lasso.R` | TCGA-COAD download (UCSC Xena), Cox, LASSO-Cox, time-dependent ROC |
| `R/04_gse39582_validation.R` | External validation in GSE39582 |
| `R/05_gse39582_subgroups.R` | Molecular subgroup analyses (original submission) |
| `R/06_supplementary_analyses.R` | VIF, xCell and Cohort-2 genes (original submission; xCell sign error fixed) |
| `R/08_revision_model_diagnostics.R` | EPV, λ.min/λ.1se, CV curve, Schoenfeld tests, collinearity, bootstrap and 100-seed stability, gene-removal sensitivity |
| `R/09_revision_platform.R` | Cross-platform harmonization (z-score, rank-INT) and reverse-direction validation |
| `R/10_revision_subgroups_utility.R` | Interaction tests, Table 3, 5-year OS with 95% CI, C-index/AUC/DeLong/NRI/IDI, decision curves, adjuvant chemotherapy |
| `R/11_revision_signatures.R` | Comparison with CMS, CIT and an Oncotype DX Colon-like score; Cohort-2 cell-cycle and curation-free analyses |
| `R/12_revision_xcell.R` | xCell in both cohorts (signed correlations), gene source, immune/stroma-adjusted Cox |
| `R/13_revision_ppi_batch.R` | PPI network statistics, top-100 hub tables, PCA before/after ComBat |
| `tools/export_processed_data.R` | Author-side script that exported `data/processed/` from the original workspace |
| `tools/add_tissue_and_cms.R` | Adds the `tissue` column to the GSE39582 table and recomputes CMS on primary tumours only (v1.0.1) |
| `data/processed/` | Derived analysis tables (see below) |
| `results/revision/` | Supplementary tables and figures produced by scripts 08–13 |

## Data

All data are public and de-identified.

- **GEO:** GSE41328, GSE37364, GSE15960 (Cohort-1), GSE21510, GSE32323, GSE110224 (Cohort-2), GSE39582 (validation). All are on platform GPL570.
- **TCGA-COAD:** GDC-harmonized STAR-FPKM and clinical data, via UCSC Xena (`UCSCXenaTools`).
- **STRING v12.0:** Homo sapiens, via `STRINGdb`, with combined score ≥ 400 and the full network.

`data/processed/` holds the derived tables used in the paper:

- TCGA (n = 436) and GSE39582 (n = 562 primary tumours + 17 non-tumoural mucosa samples, flagged in the `tissue` column and excluded by `R/00_helpers.R`) analysis tables, with clinical variables and the genes analysed
- LASSO coefficients
- xCell scores
- CMS (CMScaller) and CIT labels
- Discovery-cohort sample sheets and DEG lists
- PPI edge lists (unique undirected edges)

UCSC Xena and GDC releases change over time, so a fresh download may not reproduce n = 436 exactly. Scripts 08–13 therefore run from these tables.

## How to run

Use R ≥ 4.5 and open `crc-tme-5gene-score.Rproj` (the working directory must be the repository root).

```r
install.packages(c("survival","glmnet","timeROC","car","survminer","ggplot2","dplyr",
                   "igraph","survIDINRI","dcurves","pROC","remotes","UCSCXenaTools"))
BiocManager::install(c("GEOquery","limma","sva","clusterProfiler","org.Hs.eg.db","STRINGdb"))
remotes::install_github("dviraran/xCell")
remotes::install_github("Lothelab/CMScaller")   # only for RECOMPUTE_CMS in script 11

## Revision analyses (from data/processed; ~10 minutes)
source("run_revision.R")
```

Scripts 01–03 download the raw data (several GB) and take longer. Script 13 needs the outputs of 01–02 to draw the PCA plots.

Random seeds are fixed (`set.seed(123)` for cross-validation, `2026` for bootstrap). For `glmnet` ≥ 5.0, Breslow ties are set explicitly (`cox.ties = "breslow"`) to match the original analysis.

## Corrections relative to the submitted manuscript

1. **PPI edge counts.** STRINGdb lists each interaction twice. The unique undirected edges are 14,914 (Cohort-1) and 11,096 (Cohort-2), not 29,828 and 22,192. Node counts and hub rankings are unchanged.
2. **xCell correlations.** These were reported as absolute values. Signed values are in `R/12_revision_xcell.R`.
3. **Table 3.** Regenerated directly from the analysis output (`R/10_revision_subgroups_utility.R`).
4. **Non-tumoural samples in GSE39582 (v1.0.1).** GSE39582 contains 566 primary tumours and 19 non-tumoural mucosa samples. The first export kept 17 mucosa samples with survival data (n = 579). From v1.0.1 all analyses use primary tumours only (n = 562, 191 deaths), and CMS labels are recomputed on tumours only.

## License

MIT (see `LICENSE`). If you use this code, please cite the article above and the Zenodo archive of this repository.

## Contact

Burhanettin Yalçınkaya, TÜBİTAK UME, Bioanalysis Laboratory, Gebze, Kocaeli, Türkiye.
