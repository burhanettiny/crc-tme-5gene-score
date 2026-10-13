## ============================================================
## 03_tcga_Cox_LASSO.R
## TCGA-COAD: GDC-harmonized STAR-FPKM expression + clinical data
## Univariate Cox → Multivariate Cox → LASSO-Cox → Time-dependent ROC
## OUTPUT: tcga_final_df_GDC.rds, tcga_lasso_sonuc.rds
## NOTE (revision): UCSC Xena renamed the GDC phenotype file after the original
## analysis (June 2026); both names are accepted below. Because GDC/Xena
## releases change, a fresh download may not reproduce n = 436 exactly; the
## exact analysis table used in the paper is data/processed/tcga_coad_analysis.csv.
## ============================================================

library(UCSCXenaTools); library(dplyr)
options(timeout = 3600)   # the expression file is ~120 MB
library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)
library(survival); library(glmnet); library(timeROC)

dir.create("results", showWarnings = FALSE)
dir.create("tcga_data_gdc", showWarnings = FALSE)

## ---- 1. Download TCGA-COAD from UCSC Xena (GDC harmonized) ----
cat("Downloading TCGA-COAD expression (GDC STAR-FPKM)...\n")

q_expr <- XenaGenerate(
  subset = XenaCohorts == "GDC TCGA Colon Cancer (COAD)") %>%
  XenaFilter(filterDatasets = "star_fpkm.tsv$")
q_clin <- XenaGenerate(
  subset = XenaCohorts == "GDC TCGA Colon Cancer (COAD)") %>%
  XenaFilter(filterDatasets = "TCGA-COAD\\.(GDC_phenotype|clinical)\\.tsv$")
q_surv <- XenaGenerate(
  subset = XenaCohorts == "GDC TCGA Colon Cancer (COAD)") %>%
  XenaFilter(filterDatasets = "TCGA-COAD.survival.tsv$")

expr_df <- as.data.frame(XenaPrepare(XenaDownload(XenaQuery(q_expr))))
clin    <- as.data.frame(XenaPrepare(XenaDownload(XenaQuery(q_clin))))
surv    <- as.data.frame(XenaPrepare(XenaDownload(XenaQuery(q_surv))))

## ---- 2. Prepare expression matrix ----
rownames(expr_df) <- expr_df[[1]]
expr_mat <- as.matrix(expr_df[, -1])
rownames(expr_mat) <- sub("\\..*$", "", rownames(expr_mat))  # remove Ensembl version

## ---- 3. Map 10 TME hub genes: SYMBOL → ENSEMBL ----
tme_genes <- c("IL6","IL1B","CXCL8","MYC","PTGS2",
               "COL1A1","IL1A","TIMP1","CXCL12","MMP3")
gene_map <- bitr(tme_genes, fromType="SYMBOL", toType="ENSEMBL",
                  OrgDb=org.Hs.eg.db)
gene_map <- gene_map %>%
  filter(ENSEMBL %in% rownames(expr_mat)) %>%
  distinct(SYMBOL, .keep_all=TRUE)
cat("Genes found in expression matrix:", nrow(gene_map), "/ 10\n")

## ---- 4. Filter primary tumor samples + merge clinical ----
primary_samples <- colnames(expr_mat)[grepl("-01[A-Z]?$", colnames(expr_mat))]
cat("Primary tumor samples:", length(primary_samples), "\n")

clin_primary <- clin[clin$sample %in% primary_samples, ]
surv_primary <- surv[surv$sample %in% primary_samples, ]

stage_col  <- "ajcc_pathologic_stage.diagnoses"
age_col    <- "age_at_index.demographic"
gender_col <- "gender.demographic"

analiz_df <- clin_primary %>%
  dplyr::select(sample, stage_raw = !!stage_col,
                age   = !!age_col, gender = !!gender_col) %>%
  dplyr::inner_join(surv_primary %>% dplyr::select(sample, OS, OS.time), by="sample")

analiz_df$stage_binary <- dplyr::case_when(
  grepl("Stage I$|Stage IA|Stage IB|Stage II", analiz_df$stage_raw, ignore.case=TRUE) &
    !grepl("Stage III|Stage IV", analiz_df$stage_raw, ignore.case=TRUE) ~ "I-II",
  grepl("Stage III|Stage IV", analiz_df$stage_raw, ignore.case=TRUE) ~ "III-IV",
  TRUE ~ NA_character_)

analiz_df_complete <- analiz_df %>%
  dplyr::filter(!is.na(OS), !is.na(OS.time), !is.na(age),
                !is.na(gender), !is.na(stage_binary))
cat("Complete cases (n):", nrow(analiz_df_complete), "  [expected: 436]\n")

## ---- 5. Add gene expression ----
gene_ids <- setNames(gene_map$ENSEMBL, gene_map$SYMBOL)
expr_sub <- t(expr_mat[gene_ids, analiz_df_complete$sample, drop=FALSE])
colnames(expr_sub) <- names(gene_ids)
final_df <- dplyr::inner_join(analiz_df_complete,
                               as.data.frame(expr_sub) %>%
                                 tibble::rownames_to_column("sample"), by="sample")
cat("Final analysis table:", nrow(final_df), "x", ncol(final_df), "\n")

saveRDS(list(final_df=final_df, mevcut_genler=names(gene_ids)),
        "results/tcga_final_df_GDC.rds")

## ---- 6. Univariate Cox (High vs Low, median split) ----
cat("\n--- Univariate Cox ---\n")
for (g in names(gene_ids)) {
  med <- median(final_df[[g]], na.rm=TRUE)
  grp <- factor(ifelse(final_df[[g]] > med, "High","Low"), levels=c("Low","High"))
  cx  <- coxph(Surv(final_df$OS.time, final_df$OS) ~ grp)
  s   <- summary(cx)
  cat(sprintf("%-8s HR=%.3f p=%.4f\n", g,
              s$coefficients[1,2], s$coefficients[1,5]))
}
## Confirmed: CXCL8 HR=0.511 p=0.0025, MMP3 HR=0.573 p=0.0115,
##            PTGS2 HR=0.707 p=0.1051 (NS)

## ---- 7. Multivariate Cox (age + sex + stage) ----
cat("\n--- Multivariate Cox (age + sex + stage) ---\n")
for (g in names(gene_ids)) {
  med <- median(final_df[[g]], na.rm=TRUE)
  tmp <- final_df
  tmp$grp <- factor(ifelse(tmp[[g]] > med, "High","Low"), levels=c("Low","High"))
  tmp$stage_binary <- factor(tmp$stage_binary, levels=c("I-II","III-IV"))
  tmp$gender <- factor(tmp$gender)
  cx <- coxph(Surv(OS.time, OS) ~ grp + age + gender + stage_binary, data=tmp)
  s  <- summary(cx)
  cat(sprintf("%-8s HR=%.3f p=%.4f\n", g,
              s$coefficients["grpHigh","exp(coef)"],
              s$coefficients["grpHigh","Pr(>|z|)"]))
}
## Confirmed: CXCL8 HR=0.52 p=0.0034, MMP3 HR=0.55 p=0.0074

## ---- 8. LASSO-Cox (glmnet, 10-fold CV) ----
X <- as.matrix(final_df[, names(gene_ids)])
y <- Surv(final_df$OS.time, final_df$OS)
set.seed(123)
cvfit    <- cv.glmnet(X, y, family="cox", alpha=1, nfolds=10)
coef_min <- coef(cvfit, s="lambda.min")
cat(sprintf("\nLambda.min = %.5f  [expected: 0.00807]\n", cvfit$lambda.min))
cat("Selected genes and coefficients:\n")
print(coef_min)
## Confirmed: IL6(+0.356), TIMP1(+0.250), CXCL8(-0.066),
##            PTGS2(-0.113), MMP3(-0.162); others zeroed out

final_df$risk_score <- as.numeric(X %*% as.numeric(coef_min))

## ---- 9. KM (median split) ----
library(survminer)
final_df$risk_grp <- factor(ifelse(final_df$risk_score > median(final_df$risk_score),
                                    "High","Low"), levels=c("Low","High"))
fit_km  <- survfit(Surv(OS.time, OS) ~ risk_grp, data=final_df)
sdiff   <- survdiff(Surv(OS.time, OS) ~ risk_grp, data=final_df)
pval_km <- 1 - pchisq(sdiff$chisq, 1)
cat(sprintf("\nKM log-rank p = %.5f  [expected: 0.00029]\n", pval_km))

## ---- 10. Time-dependent ROC (3 models, 1/3/5 year) ----
times_days <- c(365, 1095, 1825)
stage_num  <- ifelse(final_df$stage_binary=="III-IV", 1, 0)
roc_stage  <- timeROC(T=final_df$OS.time, delta=final_df$OS, marker=stage_num,
                       cause=1, times=times_days, iid=FALSE)
roc_gene   <- timeROC(T=final_df$OS.time, delta=final_df$OS, marker=final_df$risk_score,
                       cause=1, times=times_days, iid=FALSE)
cx_comb    <- coxph(Surv(OS.time, OS) ~ risk_score + stage_binary, data=final_df)
final_df$combined_lp <- predict(cx_comb, type="lp")
roc_combined <- timeROC(T=final_df$OS.time, delta=final_df$OS,
                         marker=final_df$combined_lp,
                         cause=1, times=times_days, iid=FALSE)

cat("\n3-Year AUC:\n")
cat(sprintf("  Stage only:  %.3f  [expected: 0.716]\n", roc_stage$AUC[2]))
cat(sprintf("  5-gene score: %.3f  [expected: 0.628]\n", roc_gene$AUC[2]))
cat(sprintf("  Combined:    %.3f  [expected: 0.793]\n", roc_combined$AUC[2]))

saveRDS(list(final_df=final_df, cvfit=cvfit, coef_min=coef_min,
             roc_stage=roc_stage, roc_gene=roc_gene,
             roc_combined=roc_combined),
        "results/tcga_lasso_sonuc.rds")
cat("\nSaved: results/tcga_lasso_sonuc.rds\n")
