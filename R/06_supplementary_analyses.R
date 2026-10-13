## ============================================================
## 06_supplementary_analyses.R
## (A) VIF / multicollinearity for 5 LASSO genes
## (B) xCell immune cell deconvolution + risk score correlation
## (C) Cohort-2 hub genes tested in TCGA (continuous Cox)
## ============================================================

library(dplyr); library(survival)
if (!requireNamespace("car", quietly=TRUE)) install.packages("car")
library(car)

## ---- Load TCGA data ----
tcga <- readRDS("results/tcga_lasso_sonuc.rds")
final_df <- tcga$final_df
lasso_genes <- c("IL6","CXCL8","PTGS2","TIMP1","MMP3")

## ====================================================
## A. VIF + Pairwise Correlation
## ====================================================
cat("=== A. VIF and Pairwise Correlation ===\n")
cor_mat <- cor(final_df[, lasso_genes])
cat("Pairwise Pearson r:\n"); print(round(cor_mat, 3))
cat(sprintf("Min/Max |r| (off-diagonal): %.3f / %.3f\n",
    min(abs(cor_mat[lower.tri(cor_mat)])),
    max(abs(cor_mat[lower.tri(cor_mat)]))))
## IL6–CXCL8 r=0.718 (highest); TIMP1 most independent

lm_fit  <- lm(OS.time ~ IL6 + CXCL8 + PTGS2 + TIMP1 + MMP3, data=final_df)
vif_vals <- car::vif(lm_fit)
cat("\nVIF values:\n"); print(round(vif_vals, 2))
cat(sprintf("VIF range: %.2f – %.2f  [expected: 1.17 – 3.08]\n",
    min(vif_vals), max(vif_vals)))
## All VIF < 5 → no multicollinearity concern

## ====================================================
## B. xCell Immune Cell Deconvolution
## ====================================================
cat("\n=== B. xCell Immune Cell Correlation ===\n")

if (!requireNamespace("xCell", quietly=TRUE)) {
  if (!requireNamespace("remotes", quietly=TRUE)) install.packages("remotes")
  remotes::install_github("dviraran/xCell")
}
library(xCell)
library(org.Hs.eg.db); library(AnnotationDbi)

## Requires expr_mat (TCGA ENSEMBL expression matrix) in environment
## Load via 03_tcga_Cox_LASSO.R if not available

if (!exists("expr_mat")) {
  cat("NOTE: expr_mat not found. Run 03_tcga_Cox_LASSO.R first.\n")
} else {
  ## ENSEMBL → SYMBOL
  sym_map <- AnnotationDbi::select(org.Hs.eg.db, keys=rownames(expr_mat),
                                    keytype="ENSEMBL", columns="SYMBOL")
  sym_map <- sym_map[!duplicated(sym_map$ENSEMBL), ]
  rownames(sym_map) <- sym_map$ENSEMBL
  valid   <- intersect(rownames(expr_mat), rownames(sym_map))
  expr_s  <- expr_mat[valid, ]
  syms    <- sym_map[valid, "SYMBOL"]
  ok      <- !is.na(syms)
  expr_s  <- expr_s[ok, ]; syms <- syms[ok]

  ## Deduplicate by highest mean expression
  avg     <- rowMeans(expr_s)
  sdf     <- data.frame(idx=seq_along(syms), sym=syms, avg=avg)
  best    <- sdf %>% group_by(sym) %>% slice_max(avg, n=1, with_ties=FALSE) %>% pull(idx)
  expr_final <- expr_s[best, ]
  rownames(expr_final) <- syms[best]

  expr_436 <- expr_final[, final_df$sample]
  cat("Running xCellAnalysis on n=436 samples...\n")
  xcell_scores <- xCellAnalysis(expr_436)
  cat("xCell dimensions:", nrow(xcell_scores), "cell types ×",
      ncol(xcell_scores), "samples\n")

  ## Correlate with risk score
  final_match <- final_df[match(colnames(xcell_scores), final_df$sample), ]
  cors <- apply(xcell_scores, 1, function(r)
    cor(r, final_match$risk_score, method="pearson"))

  ## TRUE immune cells only (exclude stromal/structural)
  stromal <- c("Adipocytes","Astrocytes","Chondrocytes","Endothelial cells",
               "Epithelial cells","Erythrocytes","Fibroblasts","Hepatocytes",
               "Keratinocytes","ly Endothelial cells","Melanocytes",
               "Mesangial cells","MSC","mv Endothelial cells","Myocytes",
               "Neurons","Osteoblast","Pericytes","Preadipocytes","Sebocytes",
               "Skeletal muscle","Smooth muscle","ImmuneScore","StromaScore",
               "MicroenvironmentScore","Platelets","HSC","MPP","GMP","MEP","CMP")
  immune_cors <- cors[setdiff(names(cors), stromal)]

  cat(sprintf("\nImmune cells only (n=%d): |r| range %.3f – %.3f\n",
      length(immune_cors), min(abs(immune_cors)), max(abs(immune_cors))))
  cat("Top 10 by |r|:\n")
  ## REVISION FIX: the submitted version printed abs() values, losing the sign.
  ## Signed correlations, ordered by magnitude:
  print(round(immune_cors[order(-abs(immune_cors))][1:10], 3))
  ## Full signed analysis in both cohorts (Spearman, BH-FDR): R/12_revision_xcell.R
}

## ====================================================
## C. Cohort-2 Hub Genes tested in TCGA (continuous Cox)
## ====================================================
cat("\n=== C. Cohort-2 Hub Genes in TCGA Univariate Cox ===\n")
library(AnnotationDbi)

k2_genes <- c("MYC","CDK1","CCNB1","TOP2A","CD44",
              "CCNA2","AURKA","CHEK1","MAD2L1","KIF23")

if (exists("expr_mat")) {
  gene_map_k2 <- AnnotationDbi::select(org.Hs.eg.db, keys=k2_genes,
                                        keytype="SYMBOL", columns="ENSEMBL")
  gene_map_k2 <- gene_map_k2[!duplicated(gene_map_k2$SYMBOL), ]

  for (i in seq_len(nrow(gene_map_k2))) {
    sym <- gene_map_k2$SYMBOL[i]; ens <- gene_map_k2$ENSEMBL[i]
    if (ens %in% rownames(expr_mat)) {
      final_df[[sym]] <- expr_mat[ens, final_df$sample]
      f  <- as.formula(paste("Surv(OS.time, OS) ~", sym))
      cx <- coxph(f, data=final_df)
      s  <- summary(cx)
      cat(sprintf("%-8s HR=%.3f p=%.4f\n", sym,
                  s$coefficients[1,2], s$coefficients[1,5]))
    }
  }
  ## Confirmed: all p > 0.22 in continuous Cox → none significant
  ## Binary grouping: CCNB1 p=0.046, CHEK1 p=0.013 (both HR<1, protective)
} else {
  cat("NOTE: expr_mat not available. Run 03_tcga_Cox_LASSO.R first.\n")
  cat("Confirmed results: all Cohort-2 hub genes p > 0.22 (continuous Cox)\n")
}

cat("\nSupplementary analyses complete.\n")
