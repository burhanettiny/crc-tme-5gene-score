## ============================================================
## 05_gse39582_subgroups.R
## GSE39582: Mann-Whitney U, BRAF subgroup Cox/KM,
## robustness check (all molecular subgroups), FDR correction
## OUTPUT: gse39582_altgrup_sonuc.rds
## ============================================================

library(dplyr); library(survival)

## ---- 0. Load results ----
sonuc        <- readRDS("results/gse39582_genel_sonuc.rds")
clin_complete <- sonuc$clin_complete
cat("n =", nrow(clin_complete), "\n")

## ---- 1. Mann-Whitney U (5 molecular markers) ----
cat("\n--- Mann-Whitney U: Risk score by molecular status ---\n")
mw_test <- function(df, var, g1, g2) {
  d  <- df %>% dplyr::filter(!is.na(.data[[var]]))
  x1 <- d$risk_score[d[[var]] == g1]
  x2 <- d$risk_score[d[[var]] == g2]
  p  <- wilcox.test(x1, x2)$p.value
  cat(sprintf("%-10s %s(n=%d,med=%.3f) vs %s(n=%d,med=%.3f): p=%.4f\n",
      var, g1, length(x1), median(x1), g2, length(x2), median(x2), p))
  p
}
p_kras <- mw_test(clin_complete, "kras_raw", "M", "WT")
p_braf <- mw_test(clin_complete, "braf_raw", "M", "WT")  # p=0.061
p_tp53 <- mw_test(clin_complete, "tp53_raw", "M", "WT")
p_mmr  <- mw_test(clin_complete, "mmr_raw", "dMMR", "pMMR")
p_cimp <- mw_test(clin_complete, "cimp_raw", "+", "-")
cat(sprintf("Bonferroni-corrected threshold: p < %.3f\n", 0.05/5))
## Confirmed: BRAF p=0.061 (does not reach nominal threshold)

## ---- 2. BRAF subgroup analysis ----
cat("\n--- BRAF Subgroup Analysis ---\n")
for (braf_grp in c("WT","M")) {
  d <- clin_complete %>% dplyr::filter(braf_raw == braf_grp)
  cat(sprintf("\nBRAF %s (n=%d):\n", braf_grp, nrow(d)))

  ## Univariate Cox
  cx <- coxph(Surv(os_time, os_event) ~ risk_score, data=d)
  s  <- summary(cx)
  cat(sprintf("  Univariate Cox: HR=%.3f p=%.4f\n",
      s$coefficients[1,2], s$coefficients[1,5]))

  ## Tertile KM
  tq  <- quantile(d$risk_score, probs=c(1/3,2/3))
  d$rt <- cut(d$risk_score, breaks=c(-Inf,tq[1],tq[2],Inf),
               labels=c("Low","Intermediate","High"))
  sd  <- survdiff(Surv(os_time, os_event) ~ rt, data=d)
  cat(sprintf("  Log-rank p=%.4f\n", 1-pchisq(sd$chisq, 2)))

  ## Stage-adjusted
  d_s <- d %>% dplyr::filter(!is.na(stage_binary)) %>%
    dplyr::mutate(stage_binary=factor(stage_binary, levels=c("I-II","III-IV")))
  if (nrow(d_s) > 10) {
    cx2 <- coxph(Surv(os_time, os_event) ~ risk_score + stage_binary, data=d_s)
    s2  <- summary(cx2)
    cat(sprintf("  Stage-adjusted: HR=%.3f p=%.4f (n=%d)\n",
        s2$coefficients["risk_score","exp(coef)"],
        s2$coefficients["risk_score","Pr(>|z|)"],
        nrow(d_s)))
  }
}
## Confirmed: BRAF-WT HR=1.742 p=0.0016, LR p=0.0036, stage-adj p=0.0047
## Confirmed: BRAF-Mut HR=2.155 p=0.195

## ---- 3. Robustness check (Table 6) ----
cat("\n--- Robustness Check: All Molecular Subgroups ---\n")
subgroups <- list(
  list(var="kras_raw", label="KRAS-Mutant",    val="M"),
  list(var="kras_raw", label="KRAS-WT",         val="WT"),
  list(var="tp53_raw", label="TP53-Mutant",     val="M"),
  list(var="tp53_raw", label="TP53-WT",         val="WT"),
  list(var="mmr_raw",  label="MMR-pMMR",        val="pMMR"),
  list(var="mmr_raw",  label="MMR-dMMR",        val="dMMR"),
  list(var="cimp_raw", label="CIMP-Negative",   val="-"),
  list(var="cimp_raw", label="CIMP-Positive",   val="+"),
  list(var="braf_raw", label="BRAF-WT",         val="WT"),
  list(var="braf_raw", label="BRAF-Mutant",     val="M")
)

results <- data.frame()
for (sg in subgroups) {
  d <- clin_complete %>% dplyr::filter(.data[[sg$var]] == sg$val)
  if (nrow(d) < 10) next
  cx <- tryCatch(
    coxph(Surv(os_time, os_event) ~ risk_score, data=d),
    error=function(e) NULL)
  if (is.null(cx)) next
  s <- summary(cx)
  results <- rbind(results, data.frame(
    Subgroup=sg$label, n=nrow(d),
    HR=s$coefficients[1,2], p_nominal=s$coefficients[1,5]))
}

## FDR correction (Benjamini-Hochberg)
results$p_FDR_BH      <- p.adjust(results$p_nominal, method="BH")
results$p_Bonferroni  <- p.adjust(results$p_nominal, method="bonferroni")
cat("\nTable 6 with FDR correction:\n")
print(results, digits=3)
cat("\nSubgroups surviving FDR < 0.05:\n")
print(results[results$p_FDR_BH < 0.05, c("Subgroup","HR","p_nominal","p_FDR_BH")])
## Confirmed: KRAS-Mutant, MMR-pMMR, CIMP-Negative, CIMP-Positive, BRAF-WT

## ---- 4. Actual KM survival rates (Table footnote) ----
cat("\n--- 5-year survival rates (BRAF-WT subgroup) ---\n")
cc_wt <- clin_complete %>% dplyr::filter(braf_raw=="WT")
tq    <- quantile(cc_wt$risk_score, probs=c(1/3,2/3))
cc_wt$rt <- cut(cc_wt$risk_score, breaks=c(-Inf,tq[1],tq[2],Inf),
                 labels=c("Low","Intermediate","High"))
fit  <- survfit(Surv(os_time, os_event) ~ rt, data=cc_wt)
s60  <- summary(fit, times=60)
print(data.frame(Tertile=s60$strata, Surv_5yr=round(s60$surv,3),
                 CI_low=round(s60$lower,3), CI_up=round(s60$upper,3)))
## Confirmed: Low=0.695, Intermediate=0.752, High=0.581

saveRDS(list(clin_complete=clin_complete, results_table6=results,
             mw_pvals=c(KRAS=p_kras,BRAF=p_braf,TP53=p_tp53,
                        MMR=p_mmr,CIMP=p_cimp),
             roc_39582=sonuc$roc_39582),
        "results/gse39582_altgrup_sonuc.rds")
cat("\nSaved: results/gse39582_altgrup_sonuc.rds\n")
