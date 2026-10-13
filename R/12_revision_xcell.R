## ============================================================
## 12_revision_xcell.R
## Reviewer comments 5a–b and 8: xCell deconvolution in both cohorts.
## Signed Spearman correlations with the risk score (BH-FDR across all
## cell types), cellular source of the five genes, and immune/stroma-
## adjusted Cox models.
## xCell scores are stored in data/processed (xCell, TCGA: rnaseq = TRUE on
## log2(FPKM+1) gene symbols; GSE39582: rnaseq = FALSE, one probe per gene,
## highest variance). See 06_supplementary_analyses.R to regenerate them.
## NOTE: the submitted manuscript reported absolute correlations; this
## script reports signed values.
## Output: TableS15_*
## ============================================================
source("R/00_helpers.R")
rd <- function(f) { x <- read.csv(f, check.names = FALSE); m <- as.matrix(x[, -1]); rownames(m) <- x[[1]]; m }
xc_t <- rd("data/processed/xcell_tcga_coad.csv")[, tc$sample]
xc_g <- rd("data/processed/xcell_gse39582.csv")[, gs$sample]

corr <- function(xc, s) { r <- apply(xc, 1, function(v) cor.test(v, s, method = "spearman", exact = FALSE))
  data.frame(cell_type = rownames(xc), rho = sapply(r, `[[`, "estimate"), p = sapply(r, `[[`, "p.value")) }
xk <- merge(corr(xc_t, tc$risk_score), corr(xc_g, gs$risk_score), by = "cell_type", suffixes = c("_TCGA","_GSE39582"))
xk$FDR_TCGA <- p.adjust(xk$p_TCGA, "BH"); xk$FDR_GSE39582 <- p.adjust(xk$p_GSE39582, "BH")
xk$consistent_FDR05 <- xk$FDR_TCGA < 0.05 & xk$FDR_GSE39582 < 0.05 & sign(xk$rho_TCGA) == sign(xk$rho_GSE39582)
xk <- xk[order(-abs(xk$rho_TCGA)), ]
num <- sapply(xk, is.numeric); xk[num] <- lapply(xk[num], signif, 3)
print(xk[xk$consistent_FDR05, ], row.names = FALSE)
write.csv(xk, outp("TableS15_xCell_correlations.csv"), row.names = FALSE)

src <- function(xc, E) round(sapply(c("StromaScore","ImmuneScore","Epithelial cells"), function(h)
  apply(E, 2, function(x) cor(x, xc[h, ], method = "spearman"))), 2)
source_tab <- rbind(data.frame(cohort = "TCGA", gene = G5, src(xc_t, as.matrix(tc[, G5])), check.names = FALSE),
                    data.frame(cohort = "GSE39582", gene = G5, src(xc_g, as.matrix(gs[, G5])), check.names = FALSE))
print(source_tab, row.names = FALSE)
write.csv(source_tab, outp("TableS15_gene_source.csv"), row.names = FALSE)

adj <- function(xc, d, time, event, label) {
  a <- data.frame(time, event, risk = as.numeric(scale(d$risk_score)), stage = d$stage_f,
                  CD8 = xc["CD8+ T-cells", ], M2 = xc["Macrophages M2", ], Treg = xc["Tregs", ],
                  Immune = xc["ImmuneScore", ], Stroma = xc["StromaScore", ])
  a[, 5:9] <- scale(a[, 5:9])
  s <- summary(coxph(Surv(time, event) ~ risk + stage + CD8 + M2 + Treg + Immune + Stroma, data = a))
  data.frame(cohort = label, term = rownames(s$conf.int), HR = round(s$conf.int[, 1], 3),
             L95 = round(s$conf.int[, 3], 3), U95 = round(s$conf.int[, 4], 3), p = signif(s$coefficients[, 5], 3)) }
adj_tab <- rbind(adj(xc_t, tc, tc$OS.time, tc$OS, "TCGA"), adj(xc_g, gs, gs$os_time, gs$os_event, "GSE39582"))
print(adj_tab, row.names = FALSE)
write.csv(adj_tab, outp("TableS15_immune_adjusted_Cox.csv"), row.names = FALSE)
