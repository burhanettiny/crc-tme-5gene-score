## ============================================================
## 08_revision_model_diagnostics.R
## Reviewer comments 5c, 7a–c, statistical review:
## EPV, lambda.min / lambda.1se, CV curve, proportional hazards,
## stage-stratified Cox, collinearity, bootstrap and seed stability,
## gene-removal sensitivity.
## Outputs: results/revision/TableS6_*, TableS7_*, TableS14_*, FigS2, FigS3
## ============================================================
source("R/00_helpers.R")
suppressPackageStartupMessages(library(car))

## ---- 1. Sample size and EPV ----
cat("TCGA n =", nrow(tc), "| events =", sum(tc$OS),
    "| EPV (10 candidates) =", round(sum(tc$OS) / 10, 1),
    "| EPV (5 selected) =", round(sum(tc$OS) / 5, 1), "\n")
cat("GSE39582 n =", nrow(gs), "| events =", sum(gs$os_event), "\n")

## ---- 2. LASSO refit (seed 123): lambda.min vs lambda.1se ----
L <- lasso_fit(X10, y_tcga)
b1se <- as.numeric(coef(L$cv, s = "lambda.1se"))
print(data.frame(gene = TME10, original = round(b_orig[TME10], 4),
                 refit_lambda_min = round(L$b, 4), lambda_1se = round(b1se, 4)), row.names = FALSE)
cat("lambda.min =", signif(L$cv$lambda.min, 3), "(", sum(L$b != 0), "genes ) | lambda.1se =",
    signif(L$cv$lambda.1se, 3), "(", sum(b1se != 0), "genes )\n")
png(outp("FigS2_LASSO_CV_path.png"), width = 2400, height = 1100, res = 300)
par(mfrow = c(1, 2))
plot(L$cv); title("A", adj = 0, line = 2.5)
plot(L$cv$glmnet.fit, xvar = "lambda", label = TRUE)
abline(v = log(c(L$cv$lambda.min, L$cv$lambda.1se)), lty = 2); title("B", adj = 0, line = 2.5)
dev.off()

## ---- 3. Proportional hazards (Schoenfeld) and stage-stratified Cox ----
tc$gender <- factor(tc$gender)
for (g in c("CXCL8","MMP3")) tc[[paste0(g, "_grp")]] <-
  factor(ifelse(tc[[g]] > median(tc[[g]]), "High", "Low"), levels = c("Low","High"))
ph_models <- list(
  TCGA_risk_uni   = coxph(Surv(OS.time, OS) ~ risk_score, data = tc),
  TCGA_risk_multi = coxph(Surv(OS.time, OS) ~ risk_score + age + gender + stage_f, data = tc),
  TCGA_CXCL8      = coxph(Surv(OS.time, OS) ~ CXCL8_grp + age + gender + stage_f, data = tc),
  TCGA_MMP3       = coxph(Surv(OS.time, OS) ~ MMP3_grp + age + gender + stage_f, data = tc),
  GSE_risk_uni    = coxph(Surv(os_time, os_event) ~ risk_score, data = gs),
  GSE_risk_stage  = coxph(Surv(os_time, os_event) ~ risk_score + stage_f, data = gs))
ph <- do.call(rbind, lapply(names(ph_models), function(m) {
  z <- cox.zph(ph_models[[m]])$table
  data.frame(model = m, term = rownames(z), chisq = round(z[, 1], 3), p = signif(z[, 3], 3)) }))
print(ph, row.names = FALSE)
strat <- rbind(
  TCGA = summary(coxph(Surv(OS.time, OS) ~ risk_score + age + gender + strata(stage_f), data = tc))$conf.int["risk_score", ],
  GSE  = summary(coxph(Surv(os_time, os_event) ~ risk_score + strata(stage_f), data = gs))$conf.int["risk_score", ])
print(round(strat, 3))
write.csv(ph, outp("TableS7_Schoenfeld.csv"), row.names = FALSE)
write.csv(data.frame(cohort = rownames(strat), round(strat, 3)), outp("TableS7_stage_stratified_Cox.csv"), row.names = FALSE)
png(outp("FigS3_Schoenfeld.png"), width = 2400, height = 1100, res = 300)
par(mfrow = c(1, 2))
plot(cox.zph(ph_models$TCGA_risk_uni)); title("TCGA-COAD")
plot(cox.zph(ph_models$GSE_risk_uni));  title("GSE39582")
dev.off()

## ---- 4. Collinearity ----
X5 <- as.matrix(tc[, G5])
vif_v <- car::vif(lm(OS.time ~ ., data = data.frame(OS.time = tc$OS.time, X5)))
ev <- eigen(cor(X5))$values
cat("\nVIF:", round(vif_v, 2), "| condition indices:", round(sqrt(max(ev) / ev), 2), "\n")
print(round(cor(X5), 2))
unpen <- summary(coxph(Surv(OS.time, OS) ~ ., data = data.frame(OS.time = tc$OS.time, OS = tc$OS, scale(X5))))$conf.int[, c(1, 3, 4)]
print(round(unpen, 3))

## ---- 5. Bootstrap selection and sign stability (B = 500) ----
set.seed(2026); B <- 500
boot_b <- matrix(NA, B, length(TME10), dimnames = list(NULL, TME10))
for (b in 1:B) {
  i <- sample(nrow(X10), replace = TRUE)
  cv <- suppressWarnings(cvg(X10[i, ], y_tcga[i]))
  boot_b[b, ] <- as.numeric(coef(cv, s = "lambda.min"))
  if (b %% 100 == 0) cat("bootstrap", b, "\n")
}
stab <- data.frame(gene = TME10, original = round(b_orig[TME10], 3),
                   selected_pct = round(100 * colMeans(boot_b != 0), 1),
                   positive_pct = round(100 * colMeans(boot_b > 0), 1),
                   negative_pct = round(100 * colMeans(boot_b < 0), 1),
                   median_coef  = round(apply(boot_b, 2, median), 3))
print(stab, row.names = FALSE)
write.csv(stab, outp("TableS6_bootstrap_stability.csv"), row.names = FALSE)

## ---- 6. Repeated cross-validation: 100 seeds ----
seed_test <- function(X, label, S = 100) {
  res <- t(sapply(1:S, function(s) { set.seed(s); cv <- suppressWarnings(cvg(X, y_tcga))
                                     as.numeric(coef(cv, s = "lambda.min")) }))
  colnames(res) <- colnames(X); n_sel <- rowSums(res != 0)
  cat(sprintf("\n== %s == null models: %d/%d | median genes selected: %d\n", label, sum(n_sel == 0), S, median(n_sel)))
  fr <- sort(round(100 * colMeans(res != 0)), decreasing = TRUE); print(fr[fr > 0]); res }
rs10 <- seed_test(X10, "Curated 10 TME genes")
same5 <- apply(rs10 != 0, 1, function(r) setequal(colnames(rs10)[r], G5))
cat("Identical 5-gene set selected in", sum(same5), "/ 100 seeds\n")
top25 <- read.csv("data/processed/cohort1_top25_hubs.csv")$gene
cat("Top-25 Cohort-1 hub genes:", length(top25), "\n")
X25 <- as.matrix(tc[, top25])
rs25  <- seed_test(X25, "Top-25 hub genes, raw")
rs25z <- seed_test(zs(X25), "Top-25 hub genes, z-score")
write.csv(data.frame(set = c("curated10","top25_raw","top25_z"),
                     null_models = c(sum(rowSums(rs10 != 0) == 0), sum(rowSums(rs25 != 0) == 0), sum(rowSums(rs25z != 0) == 0)),
                     identical_5gene = c(sum(same5), NA, NA)),
          outp("TableS6_seed_stability.csv"), row.names = FALSE)

## ---- 7. Gene-removal sensitivity (IL6 or CXCL8 removed) ----
drop_tab <- do.call(rbind, lapply(c("IL6","CXCL8"), function(drop) {
  keep <- setdiff(TME10, drop); Ld <- lasso_fit(X10[, keep], y_tcga)
  cat("\nWithout", drop, ":", paste(names(Ld$b)[Ld$b != 0], collapse = ", "), "\n")
  rbind(evaluate_score(X10[, keep] %*% Ld$b, tc$OS.time, tc$OS, tc$stage_binary, TT_TCGA, paste("TCGA, without", drop)),
        evaluate_score(G10[, keep] %*% Ld$b, gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, paste("GSE39582, without", drop))) }))
print(drop_tab, digits = 3, row.names = FALSE)
write.csv(drop_tab, outp("TableS14_gene_removal.csv"), row.names = FALSE)
write.csv(data.frame(gene = G5, VIF = round(vif_v, 2)), outp("TableS14_VIF.csv"), row.names = FALSE)
