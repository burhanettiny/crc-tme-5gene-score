## ============================================================
## 02_cohort2_DEG_PPI.R
## Discovery Cohort-2: GSE21510 + GSE32323 + GSE110224 (GPL570), n = 94
## (36 tumour, 58 normal; sample selection in data/processed/cohort2_samples.csv)
## Duplicate patients: GSE21510 and GSE32323 share 17 patients, identified
## by matching patient identifiers in the GEO sample titles
## ("patient 006, cancer, LCM" vs "patient 006b, normal, homogenized").
## GSE32323 (paired tumour/normal) was retained for these patients.
## Same pipeline as Cohort-1. Outputs: results/cohort2_expr.rds
## ============================================================
options(timeout = 3600)
suppressPackageStartupMessages({
  library(GEOquery); library(limma); library(sva); library(STRINGdb); library(igraph)
})
dir.create("results", showWarnings = FALSE); dir.create("data/raw", recursive = TRUE, showWarnings = FALSE)

ss <- read.csv("data/processed/cohort2_samples.csv")
ss$group <- factor(ifelse(grepl("normal", ss$group, ignore.case = TRUE), "Normal", "Tumor"), levels = c("Normal","Tumor"))
print(table(ss$dataset, ss$group))

lg <- function(m) if (max(m, na.rm = TRUE) > 50) log2(pmax(m, 1)) else m
gl <- lapply(setNames(unique(ss$dataset), unique(ss$dataset)),
             function(id) getGEO(id, GSEMatrix = TRUE, destdir = "data/raw")[[1]])

## ---- 1. Documenting the 17 overlapping patients ----
pid <- function(g) unique(sub("^patient ([0-9]+).*$", "\\1", pData(g)$title[grepl("^patient [0-9]+", pData(g)$title)]))
shared <- intersect(pid(gl[["GSE21510"]]), sub("b$", "", pid(gl[["GSE32323"]])))
cat("Patients shared by GSE21510 and GSE32323:", length(shared), "\n"); print(shared)

## ---- 2. Expression, normalization, ComBat ----
ex <- lapply(gl, function(g) lg(exprs(g)))
common <- Reduce(intersect, lapply(ex, rownames))
E <- do.call(cbind, lapply(names(ex), function(i) ex[[i]][common, ss$sample[ss$dataset == i], drop = FALSE]))[, ss$sample]
E_norm   <- normalizeBetweenArrays(E, method = "quantile")
E_combat <- ComBat(E_norm, batch = factor(ss$dataset), mod = model.matrix(~ ss$group))
saveRDS(list(norm = E_norm, combat = E_combat, batch = factor(ss$dataset), group = ss$group),
        "results/cohort2_expr.rds")

## ---- 3. Differential expression ----
d   <- model.matrix(~ 0 + ss$group); colnames(d) <- levels(ss$group)
fit <- eBayes(contrasts.fit(lmFit(E_combat, d), makeContrasts(Tumor - Normal, levels = d)))
tt  <- topTable(fit, number = Inf, adjust.method = "BH")
sig <- tt[tt$adj.P.Val < 0.05 & abs(tt$logFC) > 1, ]
orig <- read.csv("data/processed/cohort2_DEG.csv")
cat("DEGs (re-run):", nrow(sig), "| published list:", nrow(orig),
    "| published probes recovered:", sum(orig$probe %in% rownames(sig)), "\n")

## ---- 4. PPI network (published DEG list) ----
sym_col <- intersect(c("gene_symbol","gene","SYMBOL"), names(orig))[1]
genes <- unique(na.omit(orig[[sym_col]])); genes <- genes[genes != ""]
cat("Unique gene symbols:", length(genes), "(published: 1,337)\n")
sdb <- STRINGdb$new(version = "12.0", species = 9606, score_threshold = 400,
                    network_type = "full", input_directory = "data/raw")
mp  <- sdb$map(data.frame(gene = genes), "gene", removeUnmappedRows = TRUE)
int <- sdb$get_interactions(unique(mp$STRING_id))
g   <- simplify(graph_from_data_frame(int[, 1:2], directed = FALSE))
cat("STRINGdb rows:", nrow(int), "| nodes:", vcount(g), "| unique edges:", ecount(g), "\n")
sym <- setNames(mp$gene, mp$STRING_id); sym <- sym[!duplicated(names(sym))]
print(head(sort(setNames(degree(g), sym[V(g)$name]), decreasing = TRUE), 10))
