## ============================================================
## 09_revision_platform.R
## Reviewer comment 2: cross-platform harmonization.
## (i) within-cohort z-score, (ii) rank-based inverse normal
## transformation (rank-INT), (iii) reverse direction
## (train on GSE39582 microarray, test on TCGA RNA-seq).
## Output: results/revision/TableS8_platform_harmonization.csv
## ============================================================
source("R/00_helpers.R")

L_z   <- lasso_fit(zs(X10),  y_tcga)
L_int <- lasso_fit(int(X10), y_tcga)
print(round(cbind(original = b_orig[TME10], z_score = L_z$b, rank_INT = L_int$b), 4))

s_z   <- zs(G10)  %*% L_z$b
s_int <- int(G10) %*% L_int$b
cat("Spearman rho with original score (GSE39582): z =", round(cor(gs$risk_score, s_z, method = "s"), 3),
    "| rank-INT =", round(cor(gs$risk_score, s_int, method = "s"), 3), "\n")

## Reverse direction: patients with zero follow-up cannot enter glmnet Cox (training only)
pos   <- gs$os_time > 0
cat("GSE39582 patients with zero follow-up (excluded from training):", sum(!pos), "\n")
L_gse <- lasso_fit(zs(G10)[pos, ], Surv(gs$os_time[pos], gs$os_event[pos]))
cat("GSE39582-trained coefficients:\n"); print(round(L_gse$b[L_gse$b != 0], 4))
cat("lambda.1se genes:", sum(coef(L_gse$cv, s = "lambda.1se") != 0), "\n")

tab <- rbind(
  evaluate_score(gs$risk_score, gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "GSE39582: original (unharmonized)"),
  evaluate_score(s_z,   gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "GSE39582: z-score harmonized"),
  evaluate_score(s_int, gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "GSE39582: rank-INT harmonized"),
  evaluate_score(tc$risk_score, tc$OS.time, tc$OS, tc$stage_binary, TT_TCGA, "TCGA: original (training, apparent)"),
  evaluate_score(zs(X10) %*% L_gse$b, tc$OS.time, tc$OS, tc$stage_binary, TT_TCGA, "TCGA: GSE39582-trained score (reverse, z)"))
num <- sapply(tab, is.numeric); tab[num] <- lapply(tab[num], signif, 3)
print(tab, row.names = FALSE)
write.csv(tab, outp("TableS8_platform_harmonization.csv"), row.names = FALSE)
write.csv(data.frame(gene = TME10, original = b_orig[TME10], z_score = L_z$b, rank_INT = L_int$b,
                     GSE_trained_z = L_gse$b), outp("TableS8_coefficients.csv"), row.names = FALSE)
