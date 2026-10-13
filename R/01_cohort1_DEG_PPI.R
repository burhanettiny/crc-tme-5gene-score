## ============================================================
## 01_cohort1_DEG_PPI.R
## Discovery Cohort-1: GSE41328 + GSE37364 + GSE15960 (GPL570), n = 84
## (30 tumour, 54 normal; sample selection in data/processed/cohort1_samples.csv)
## log2 -> quantile normalization -> ComBat (batch = dataset,
## mod = ~group) -> limma (|log2FC| > 1, BH-adjusted p < 0.05)
## -> GO/KEGG -> STRING v12.0 PPI (STRINGdb, score >= 400, full network)
## Outputs: results/cohort1_expr.rds, results/cohort1_DEG_rerun.csv
## ============================================================
options(timeout = 3600)
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(sva); library(clusterProfiler)
  library(org.Hs.eg.db); library(STRINGdb); library(igraph)
})
dir.create("results", showWarnings = FALSE); dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

ss <- read.csv("data/processed/cohort1_samples.csv")
ss$group <- factor(ss$group, levels = c("Normal","Tumor"))
print(table(ss$dataset, ss$group))

## ---- 1. Download and select samples ----
lg <- function(m) if (max(m, na.rm = TRUE) > 50) log2(pmax(m, 1)) else m      # log2 only if not yet log-scaled
get_expr <- function(id) { g <- getGEO(id, GSEMatrix = TRUE, destdir = "data/raw")[[1]]; lg(exprs(g)) }
ids <- unique(ss$dataset)
ex  <- lapply(setNames(ids, ids), get_expr)
common <- Reduce(intersect, lapply(ex, rownames))
E <- do.call(cbind, lapply(ids, function(i) ex[[i]][common, ss$sample[ss$dataset == i], drop = FALSE]))
E <- E[, ss$sample]

## ---- 2. Normalization and batch correction ----
E_norm   <- normalizeBetweenArrays(E, method = "quantile")
E_combat <- ComBat(E_norm, batch = factor(ss$dataset), mod = model.matrix(~ ss$group))
saveRDS(list(norm = E_norm, combat = E_combat, batch = factor(ss$dataset), group = ss$group),
        "results/cohort1_expr.rds")

## ---- 3. Differential expression ----
d  <- model.matrix(~ 0 + ss$group); colnames(d) <- levels(ss$group)
fit <- eBayes(contrasts.fit(lmFit(E_combat, d), makeContrasts(Tumor - Normal, levels = d)))
tt <- topTable(fit, number = Inf, adjust.method = "BH")
sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, ]
orig <- read.csv("data/processed/cohort1_DEG_2969.csv")
cat("DEGs (re-run):", nrow(sig), "| published list:", nrow(orig),
    "| published probes recovered:", sum(orig$probe %in% rownames(sig)), "\n")
write.csv(data.frame(probe = rownames(sig), sig), "results/cohort1_DEG_rerun.csv", row.names = FALSE)

## ---- 4. Enrichment (published DEG list) ----
sym_col <- intersect(c("gene_symbol","gene","SYMBOL"), names(orig))[1]
genes <- unique(na.omit(orig[[sym_col]])); genes <- genes[genes != ""]
cat("Unique gene symbols:", length(genes), "(published: 2,002)\n")
entrez <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
ego   <- enrichGO(entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, readable = TRUE)
ekegg <- enrichKEGG(entrez, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05)
write.csv(as.data.frame(ego),   "results/cohort1_GO_BP.csv", row.names = FALSE)
write.csv(as.data.frame(ekegg), "results/cohort1_KEGG.csv",  row.names = FALSE)

## ---- 5. PPI network: STRING v12.0 via STRINGdb ----
## network_type = "full", combined score >= 400, all evidence channels (default).
## STRINGdb returns every interaction twice; unique undirected edges are counted.
sdb <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400,
                    network_type = "full", input_directory = "data/raw")
mp  <- sdb$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
int <- sdb$get_interactions(unique(mp$STRING_id))
g   <- simplify(graph_from_data_frame(int[, 1:2], directed = FALSE))
cat("STRINGdb rows:", nrow(int), "| nodes:", vcount(g), "| unique edges:", ecount(g), "\n")
sym <- setNames(mp$gene, mp$STRING_id); sym <- sym[!duplicated(names(sym))]
dg  <- sort(setNames(degree(g), sym[V(g)$name]), decreasing = TRUE)
cat("Top 25 nodes by degree:\n"); print(head(dg, 25))
## Hub definition: top 5% by degree (here 83 nodes, degree >= 63). The 10 TME-related genes
## analysed downstream were selected by the investigator from the top-ranked nodes on biological
## grounds (cytokine/chemokine signalling and extracellular-matrix remodelling); all rank within
## the top 25. A curation-free sensitivity analysis is in R/08 and R/11.
