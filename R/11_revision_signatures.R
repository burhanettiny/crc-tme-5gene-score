## ============================================================
## 11_revision_signatures.R
## Reviewer comments 1 and 3c:
##  (A) comparison with CMS (CMScaller), CIT subtypes and a microarray
##      surrogate of the Oncotype DX Colon Recurrence Score (GSE39582);
##  (B) parallel LASSO-Cox for the Cohort-2 cell-cycle hub genes,
##      outcome-independent pathway scores, penalized-spline test;
##  (C) curation-free model built from the top-25 Cohort-1 hub genes.
## CMS labels are stored in data/processed (CMScaller v2.0.1, nearest
## template prediction, 1000 permutations, FDR<0.05; set
## RECOMPUTE_CMS <- TRUE to regenerate them from GEO).
## Outputs: TableS12_*, TableS13_*
## ============================================================
source("R/00_helpers.R")
RECOMPUTE_CMS <- FALSE
cc <- gs

## ---- A1. CMS labels ----
if (RECOMPUTE_CMS) {
  suppressPackageStartupMessages({ library(GEOquery); library(CMScaller) })
  g <- getGEO("GSE39582", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]]
  E <- exprs(g)[, cc$sample]; ent <- as.character(fData(g)[["Gene ID"]])
  ok <- !is.na(ent) & ent != "" & !grepl("///", ent); E <- E[ok, ]; ent <- ent[ok]
  vv <- apply(E, 1, var); keep <- tapply(seq_along(ent), ent, function(i) i[which.max(vv[i])])
  E <- E[unlist(keep), ]; rownames(E) <- names(keep)
  set.seed(123)
  res <- CMScaller(t(scale(t(E))), rowNames = "entrez", RNAseq = FALSE, nPerm = 1000, FDR = 0.05, doPlot = FALSE)
  lab <- as.character(res$prediction); lab[is.na(lab)] <- "Unclassified"
  cc$CMS <- lab[match(cc$sample, rownames(res))]
}
cc$CMS <- factor(ifelse(is.na(cc$CMS), "Unclassified", cc$CMS),
                 levels = c("CMS1","CMS2","CMS3","CMS4","Unclassified"))
cc$CIT <- factor(cc$CIT)
print(table(cc$CMS)); print(table(cc$CIT, useNA = "ifany"))
cat("Kruskal-Wallis (risk score ~ CMS): p =", signif(kruskal.test(risk_score ~ CMS, data = cc)$p.value, 3), "\n")
print(aggregate(risk_score ~ CMS, data = cc, median))

## ---- A2. Oncotype DX Colon-like surrogate (published genes and weights) ----
refm <- rowMeans(cc[, c("ATP5E","GPX1","PGK1","UBB","VDAC2")])
nz <- function(g) cc[[g]] - refm
stromal <- (nz("BGN") + nz("FAP") + nz("INHBA")) / 3
cellcyc <- (nz("MKI67") + nz("MYC") + nz("MYBL2")) / 3
cc$RS_onco <- 0.15 * stromal - 0.30 * cellcyc + 0.15 * nz("GADD45B")
cat("Spearman: 5-gene ~ Oncotype-like =", round(cor(cc$risk_score, cc$RS_onco, method = "s"), 3),
    "| 5-gene ~ stromal group =", round(cor(cc$risk_score, stromal, method = "s"), 3), "\n")

## ---- A3. Nested Cox models on a common complete-case set ----
d5 <- cc[!is.na(cc$CIT) & !is.na(cc$stage_f), ]
cat("Common complete-case set: n =", nrow(d5), "| events =", sum(d5$os_event), "\n")
mods <- list(
  "Stage" = ~ stage_f, "Stage + 5-gene" = ~ stage_f + risk_score,
  "Stage + Oncotype-like" = ~ stage_f + RS_onco, "Stage + CMS" = ~ stage_f + CMS, "Stage + CIT" = ~ stage_f + CIT,
  "Stage + Oncotype-like + 5-gene" = ~ stage_f + RS_onco + risk_score,
  "Stage + CMS + 5-gene" = ~ stage_f + CMS + risk_score, "Stage + CIT + 5-gene" = ~ stage_f + CIT + risk_score)
fits <- lapply(mods, function(f) coxph(update(f, Surv(os_time, os_event) ~ .), data = d5))
hr5 <- t(sapply(fits, function(f) { s <- summary(f)
  if (!"risk_score" %in% rownames(s$coefficients)) return(rep(NA, 4))
  c(s$conf.int["risk_score", c(1, 3, 4)], s$coefficients["risk_score", 5]) }))
lrt <- c("Stage + 5-gene" = "Stage", "Stage + Oncotype-like + 5-gene" = "Stage + Oncotype-like",
         "Stage + CMS + 5-gene" = "Stage + CMS", "Stage + CIT + 5-gene" = "Stage + CIT")
cmp <- data.frame(model = names(fits),
                  C_index = sapply(fits, function(f) round(concordance(f)$concordance, 3)),
                  AIC = sapply(fits, function(f) round(AIC(f), 1)),
                  HR_5gene = round(hr5[, 1], 2), L95 = round(hr5[, 2], 2), U95 = round(hr5[, 3], 2),
                  p_5gene = signif(hr5[, 4], 3),
                  LRT_p_adding_5gene = sapply(names(fits), function(m)
                    if (m %in% names(lrt)) signif(anova(fits[[lrt[m]]], fits[[m]])$`Pr(>|Chi|)`[2], 3) else NA))
print(cmp, row.names = FALSE)
write.csv(cmp, outp("TableS12_signature_comparison.csv"), row.names = FALSE)

within <- do.call(rbind, lapply(levels(cc$CMS), function(k) { d <- cc[cc$CMS == k, ]
  s <- summary(coxph(Surv(os_time, os_event) ~ risk_score, data = d))
  data.frame(CMS = k, n = nrow(d), events = sum(d$os_event), HR = s$conf.int[1, 1],
             L95 = s$conf.int[1, 3], U95 = s$conf.int[1, 4], p = s$coefficients[1, 5]) }))
p_cms_int <- anova(coxph(Surv(os_time, os_event) ~ risk_score + CMS, data = cc),
                   coxph(Surv(os_time, os_event) ~ risk_score * CMS, data = cc))$`Pr(>|Chi|)`[2]
print(within, digits = 3, row.names = FALSE); cat("Risk score x CMS interaction p =", signif(p_cms_int, 3), "\n")
write.csv(within, outp("TableS12_within_CMS.csv"), row.names = FALSE)

## ---- B. Cohort-2 cell-cycle genes ----
Xk2 <- as.matrix(tc[, K2]); Gk2 <- as.matrix(gs[, K2])
single <- function(X, time, event, label) data.frame(cohort = label, gene = colnames(X), t(sapply(colnames(X), function(g) {
  s <- summary(coxph(Surv(time, event) ~ scale(X[, g]))); c(HR_SD = s$conf.int[1, 1], p = s$coefficients[1, 5]) })))
sg <- rbind(single(Xk2, tc$OS.time, tc$OS, "TCGA"), single(Gk2, gs$os_time, gs$os_event, "GSE39582"))
print(sg, digits = 3, row.names = FALSE)
L_k2 <- lasso_fit(Xk2, y_tcga)
cat("Cohort-2 LASSO: genes selected at lambda.min =", sum(L_k2$b != 0),
    "| lambda.1se =", sum(coef(L_k2$cv, s = "lambda.1se") != 0), "\n")
mz <- function(M) rowMeans(scale(M))
path_tab <- rbind(
  evaluate_score(mz(X10), tc$OS.time, tc$OS, tc$stage_binary, TT_TCGA, "TME mean z | TCGA"),
  evaluate_score(mz(G10), gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "TME mean z | GSE39582"),
  evaluate_score(mz(Xk2), tc$OS.time, tc$OS, tc$stage_binary, TT_TCGA, "Cell-cycle mean z | TCGA"),
  evaluate_score(mz(Gk2), gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "Cell-cycle mean z | GSE39582"))
for (z in list(list(mz(Xk2), tc$OS.time, tc$OS, "TCGA"), list(mz(Gk2), gs$os_time, gs$os_event, "GSE39582"))) {
  s <- z[[1]]; f <- coxph(Surv(z[[2]], z[[3]]) ~ pspline(s, df = 4))
  cat("Cell-cycle score, penalized spline,", z[[4]], ": p non-linear =",
      signif(summary(f)$coefficients[2, "p"], 3), "\n") }

## ---- C. Curation-free top-25 Cohort-1 hub genes ----
top25 <- read.csv("data/processed/cohort1_top25_hubs.csv")$gene
X25 <- as.matrix(tc[, top25]); G25 <- as.matrix(gs[, top25])
L25z <- lasso_fit(zs(X25), y_tcga)
cat("Top-25 (z) genes selected:", paste(names(L25z$b)[L25z$b != 0], collapse = ", "), "\n")
path_tab <- rbind(path_tab,
  evaluate_score(zs(G25) %*% L25z$b, gs$os_time, gs$os_event, gs$stage_binary, TT_GSE, "Top-25 hubs, curation-free (z) | GSE39582"))
num <- sapply(path_tab, is.numeric); path_tab[num] <- lapply(path_tab[num], signif, 3)
print(path_tab, row.names = FALSE)
write.csv(sg, outp("TableS13_cellcycle_single_genes.csv"), row.names = FALSE)
write.csv(path_tab, outp("TableS13_pathway_and_curation_free.csv"), row.names = FALSE)
