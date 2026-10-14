## ============================================================
## 04_gse39582_validation.R
## GSE39582 external validation cohort (Marisa et al., 2013)
## Probe collapsing → risk score → univariate/multivariate Cox
## → tertile KM → time-dependent ROC
## OUTPUT: gse39582_genel_sonuc.rds
## ============================================================

library(GEOquery); library(dplyr); library(survival); library(timeROC)

dir.create("results", showWarnings = FALSE)
dir.create("geo_data_gse39582", showWarnings = FALSE)

## ---- 1. Download GSE39582 ----
cat("Downloading GSE39582 (Marisa et al., 2013)...\n")
g <- getGEO("GSE39582", GSEMatrix=TRUE, AnnotGPL=TRUE,
             destdir="geo_data_gse39582")
gset    <- g[[1]]
pdat    <- pData(gset)
exprs_mat <- exprs(gset)
fdat    <- fData(gset)
cat("Samples:", nrow(pdat), " | Probes:", nrow(exprs_mat), "\n")

## ---- 2. Extract clinical variables ----
to_na <- function(x) ifelse(x %in% c("N/A","NA",""), NA, x)

clin39582 <- data.frame(
  sample     = rownames(pdat),
  sex        = pdat[["Sex:ch1"]],
  age        = as.numeric(pdat[["age.at.diagnosis (year):ch1"]]),
  tnm_stage  = to_na(pdat[["tnm.stage:ch1"]]),
  os_event   = as.numeric(to_na(pdat[["os.event:ch1"]])),
  os_time    = as.numeric(to_na(pdat[["os.delay (months):ch1"]])),
  mmr_raw    = to_na(pdat[["mmr.status:ch1"]]),
  cimp_raw   = to_na(pdat[["cimp.status:ch1"]]),
  tp53_raw   = to_na(pdat[["tp53.mutation:ch1"]]),
  kras_raw   = to_na(pdat[["kras.mutation:ch1"]]),
  braf_raw   = to_na(pdat[["braf.mutation:ch1"]]),
  stringsAsFactors = FALSE)

clin39582$stage_binary <- dplyr::case_when(
  clin39582$tnm_stage %in% c("0","1","2") ~ "I-II",
  clin39582$tnm_stage %in% c("3","4")     ~ "III-IV",
  TRUE ~ NA_character_)

## Complete cases for survival
## REVISION FIX: GSE39582 also contains 19 non-tumoural mucosa samples; keep primary tumours only
clin39582$tissue <- ifelse(grepl("Adenocarcinoma", pdat[["source_name_ch1"]]), "tumor", "normal")
clin_complete <- clin39582 %>%
  dplyr::filter(tissue == "tumor", !is.na(os_event), !is.na(os_time))
cat("Complete cases, primary tumours (n):", nrow(clin_complete), "  [expected: 562]\n")

## ---- 3. Probe collapsing (highest variance per gene) ----
gene_col <- grep("Gene symbol", colnames(fdat), ignore.case=TRUE, value=TRUE)[1]
probe_map <- setNames(fdat[[gene_col]], rownames(fdat))

## LASSO coefficients from TCGA (λ=0.0081)
tcga_coef  <- c(IL6=0.35577777, CXCL8=-0.06555499, PTGS2=-0.11321592,
                TIMP1=0.24952606, MMP3=-0.16229026)
lasso_genes <- names(tcga_coef)

secilen_probe <- list()
for (g_name in lasso_genes) {
  probes_g <- names(probe_map)[probe_map == g_name]
  probes_g <- probes_g[probes_g %in% rownames(exprs_mat)]
  if (!length(probes_g)) { cat("WARNING:", g_name, "not found\n"); next }
  vars <- apply(exprs_mat[probes_g,,drop=FALSE], 1, var, na.rm=TRUE)
  secilen_probe[[g_name]] <- names(which.max(vars))
  cat(sprintf("%-6s → %d probe(s), selected: %s (var=%.3f)\n",
              g_name, length(probes_g), secilen_probe[[g_name]], max(vars)))
}

## ---- 4. Calculate risk score ----
expr_genes <- sapply(lasso_genes, function(g) exprs_mat[secilen_probe[[g]], ])
colnames(expr_genes) <- lasso_genes
rownames(expr_genes) <- colnames(exprs_mat)

risk_all <- as.numeric(expr_genes %*% tcga_coef[lasso_genes])
names(risk_all) <- rownames(expr_genes)
clin_complete$risk_score <- risk_all[clin_complete$sample]
cat("\nRisk score summary:\n"); print(summary(clin_complete$risk_score))

## ---- 5. Univariate Cox (continuous risk score) ----
cx_uni <- coxph(Surv(os_time, os_event) ~ risk_score, data=clin_complete)
s_uni  <- summary(cx_uni)
cat(sprintf("\nUnivariate Cox: HR=%.3f p=%.4f  [expected: HR=1.670 p=0.0014]\n",
    s_uni$coefficients[1,2], s_uni$coefficients[1,5]))

## ---- 6. Stage-adjusted multivariate Cox ----
clin_stage <- clin_complete %>%
  dplyr::filter(!is.na(stage_binary)) %>%
  dplyr::mutate(stage_binary = factor(stage_binary, levels=c("I-II","III-IV")))
cx_multi <- coxph(Surv(os_time, os_event) ~ risk_score + stage_binary,
                   data=clin_stage)
s_multi  <- summary(cx_multi)
cat(sprintf("Stage-adjusted: HR=%.3f p=%.4f  [expected: HR=1.577 p=0.0040]\n",
    s_multi$coefficients["risk_score","exp(coef)"],
    s_multi$coefficients["risk_score","Pr(>|z|)"]))
cat(sprintf("Stage HR=%.3f p=%.4f  [expected: HR=1.723 p=0.0002]\n",
    s_multi$coefficients["stage_binaryIII-IV","exp(coef)"],
    s_multi$coefficients["stage_binaryIII-IV","Pr(>|z|)"]))

## ---- 7. Tertile KM ----
terts <- quantile(clin_complete$risk_score, probs=c(1/3, 2/3))
clin_complete$risk_tertile <- cut(clin_complete$risk_score,
                                   breaks=c(-Inf, terts[1], terts[2], Inf),
                                   labels=c("Low","Intermediate","High"))
sdiff  <- survdiff(Surv(os_time, os_event) ~ risk_tertile, data=clin_complete)
pval_lr <- 1 - pchisq(sdiff$chisq, 2)
cat(sprintf("Tertile log-rank p=%.4f  [expected: 0.013]\n", pval_lr))

## ---- 8. Time-dependent ROC (1/3/5 year; months) ----
roc_39582 <- timeROC(T=clin_complete$os_time, delta=clin_complete$os_event,
                      marker=clin_complete$risk_score, cause=1,
                      times=c(12, 36, 60), iid=TRUE)
ci_obj <- confint(roc_39582, level=0.95)
cat(sprintf("3-year AUC=%.3f (95%% CI: %.3f–%.3f)  [expected: 0.587 (0.525–0.650)]\n",
    roc_39582$AUC[2],
    ci_obj$CI_AUC["t=36","2.5%"]/100,
    ci_obj$CI_AUC["t=36","97.5%"]/100))

## ---- 9. Individual gene Cox (continuous, Table 4) ----
cat("\n--- Individual gene Cox (continuous expression) ---\n")
for (g_name in lasso_genes) {
  clin_complete[[g_name]] <- exprs_mat[secilen_probe[[g_name]],
                                        clin_complete$sample]
  f  <- as.formula(paste("Surv(os_time, os_event) ~", g_name))
  cx <- coxph(f, data=clin_complete)
  s  <- summary(cx)
  cat(sprintf("%-6s HR=%.3f p=%.4f\n", g_name,
              s$coefficients[1,2], s$coefficients[1,5]))
}
## Confirmed: TIMP1 HR=1.225 p=0.022 (only significant individual gene)

saveRDS(list(clin_complete=clin_complete, secilen_probe=secilen_probe,
             tcga_coef=tcga_coef, cx_uni=cx_uni, cx_multi=cx_multi,
             roc_39582=roc_39582),
        "results/gse39582_genel_sonuc.rds")
cat("\nSaved: results/gse39582_genel_sonuc.rds\n")
