## ============================================================
## 13_revision_ppi_batch.R
## Reviewer comment 6: PPI network statistics, hub ranks and batch
## correction diagnostics.
## NOTE: STRINGdb lists every interaction twice (A-B and B-A). The
## submitted manuscript reported these raw row counts as edges
## (29,828 and 22,192); unique undirected edges are 14,914 and 11,096.
## PCA requires the normalized matrices written by scripts 01 and 02
## (results/cohort1_expr.rds, results/cohort2_expr.rds).
## Outputs: TableS2–S5, FigS4, FigS5
## ============================================================
suppressPackageStartupMessages(library(igraph))
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)
outp <- function(f) file.path("results/revision", f)
TME10 <- c("IL6","IL1B","CXCL8","MYC","PTGS2","COL1A1","IL1A","TIMP1","CXCL12","MMP3")

net <- function(f) {
  e <- read.csv(f)
  g <- simplify(graph_from_data_frame(e[, c("string_id1","string_id2")], directed = FALSE))
  sym <- setNames(c(e$gene1, e$gene2), c(e$string_id1, e$string_id2)); sym <- sym[!duplicated(names(sym))]
  V(g)$gene <- ifelse(is.na(sym[V(g)$name]), V(g)$name, sym[V(g)$name]); g }
g1 <- net("data/processed/cohort1_ppi_edges.csv"); g2 <- net("data/processed/cohort2_ppi_edges.csv")

summ <- function(g, label) {
  dg <- degree(g); pl <- fit_power_law(dg[dg > 0], implementation = "plfit", p.value = TRUE)
  data.frame(network = label, nodes = vcount(g), unique_edges = ecount(g), density = signif(edge_density(g), 3),
             components = igraph::components(g)$no, largest_component = max(igraph::components(g)$csize),
             mean_degree = round(mean(dg), 1), median_degree = median(dg),
             global_clustering = round(transitivity(g, "global"), 3),
             powerlaw_alpha = round(pl$alpha, 2), powerlaw_xmin = pl$xmin, powerlaw_KS_p = signif(pl$KS.p, 3),
             top5pct_degree_threshold = unname(quantile(dg, 0.95)), top5pct_nodes = sum(dg >= quantile(dg, 0.95))) }
s_tab <- rbind(summ(g1, "Cohort-1"), summ(g2, "Cohort-2")); print(t(s_tab))
write.csv(s_tab, outp("TableS2_PPI_network_statistics.csv"), row.names = FALSE)

top <- function(g, n = 100) {
  t <- data.frame(gene = V(g)$gene, degree = degree(g), betweenness = round(betweenness(g, normalized = TRUE), 4),
                  closeness = round(closeness(g, normalized = TRUE), 4),
                  clustering = round(transitivity(g, type = "local", isolates = "zero"), 3))
  t <- t[order(-t$degree, -t$betweenness), ]; t$rank <- seq_len(nrow(t)); head(t, n) }
t1 <- top(g1); t2 <- top(g2)
print(t1[t1$gene %in% TME10, ], row.names = FALSE)
write.csv(t1, outp("TableS3_top100_Cohort1.csv"), row.names = FALSE)
write.csv(t2, outp("TableS4_top100_Cohort2.csv"), row.names = FALSE)

png(outp("FigS4_PPI_degree_distribution.png"), width = 2400, height = 1100, res = 300)
par(mfrow = c(1, 2))
for (x in list(list(g1, "Cohort-1"), list(g2, "Cohort-2"))) {
  dd <- degree_distribution(x[[1]], cumulative = TRUE); k <- seq_along(dd) - 1
  plot(k[k > 0], dd[k > 0], log = "xy", pch = 16, cex = 0.6, main = x[[2]], xlab = "Degree (k)", ylab = "P(K >= k)")
  abline(v = quantile(degree(x[[1]]), 0.95), lty = 2, col = "red") }
dev.off()

## ---- PCA before / after ComBat ----
pca_row <- function(E, batch, grp, title) {
  E <- E[apply(E, 1, var) > 0, ]; p <- prcomp(t(E)); ve <- round(100 * p$sdev^2 / sum(p$sdev^2), 1)
  r2 <- function(k, f) summary(lm(p$x[, k] ~ f))$r.squared
  plot(p$x[, 1], p$x[, 2], col = as.integer(factor(batch)) + 1, pch = ifelse(grp == "Normal", 1, 17),
       xlab = paste0("PC1 (", ve[1], "%)"), ylab = paste0("PC2 (", ve[2], "%)"), main = title, cex = 0.9)
  data.frame(matrix = title, PC1_var = ve[1], PC2_var = ve[2],
             R2_batch_PC1 = round(r2(1, batch), 3), R2_batch_PC2 = round(r2(2, batch), 3),
             R2_group_PC1 = round(r2(1, grp), 3), R2_group_PC2 = round(r2(2, grp), 3)) }
if (all(file.exists(c("results/cohort1_expr.rds", "results/cohort2_expr.rds")))) {
  c1 <- readRDS("results/cohort1_expr.rds"); c2 <- readRDS("results/cohort2_expr.rds")
  png(outp("FigS5_PCA_ComBat.png"), width = 2600, height = 2400, res = 300)
  par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1))
  pca_tab <- rbind(pca_row(c1$norm,   c1$batch, c1$group, "Cohort-1: before ComBat"),
                   pca_row(c1$combat, c1$batch, c1$group, "Cohort-1: after ComBat"),
                   pca_row(c2$norm,   c2$batch, c2$group, "Cohort-2: before ComBat"),
                   pca_row(c2$combat, c2$batch, c2$group, "Cohort-2: after ComBat"))
  dev.off(); print(pca_tab, row.names = FALSE)
  write.csv(pca_tab, outp("TableS5_PCA_batch.csv"), row.names = FALSE)
} else message("Run R/01 and R/02 first to create results/cohort1_expr.rds and results/cohort2_expr.rds")
