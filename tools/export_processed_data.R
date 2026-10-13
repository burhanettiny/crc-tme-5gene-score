## ============================================================
## tools/export_processed_data.R
## Run ONCE in the author's original R workspace (.RData loaded,
## plus objects created during the revision: cc, G10, Xk2, Gk2,
## X25, G25, xc_t, xc_g, int1, mp1, gp, E1c, b1, grp1).
## Writes the small, derived analysis tables to data/processed/
## so that every revision analysis can be re-run without the
## multi-GB raw downloads.
## ============================================================

## >>> Set this to the local clone of the GitHub repository <<<
REPO <- "D:/1-Kisisel/1-Akademik/5-Yayinlar/29-Col-Can/github/crc-tme-5gene-score"
out  <- file.path(REPO, "data", "processed")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
wcsv <- function(x, f) { write.csv(x, file.path(out, f), row.names = FALSE); cat("  wrote", f, dim(x), "\n") }

## ---- 1. TCGA-COAD analysis table (n = 436) ----
tc <- final_df[, c("sample","age","gender","stage_raw","stage_binary","OS","OS.time")]
tc <- cbind(tc, as.data.frame(as.matrix(final_df[, rownames(coef_min)])))          # 10 TME genes
add_cols <- function(base, M) { M <- as.data.frame(M); M[, setdiff(colnames(M), colnames(base)), drop = FALSE] }
tc <- cbind(tc, add_cols(tc, Xk2), add_cols(tc, X25))                                # Cohort-2 + top-25 hubs
cl <- clin[match(tc$sample, clin$sample), ]
tc$race   <- cl$race.demographic
tc$m_stage <- cl$ajcc_pathologic_m.diagnoses
stopifnot(nrow(tc) == 436)
wcsv(tc, "tcga_coad_analysis.csv")
wcsv(data.frame(gene = rownames(coef_min), coef = as.numeric(coef_min)), "tcga_lasso_coefficients.csv")
wcsv(data.frame(rank = 1:ncol(X25), gene = colnames(X25)), "cohort1_top25_hubs.csv")

## ---- 2. GSE39582 analysis table (n = 579) ----
onco <- c("BGN","FAP","INHBA","MKI67","MYC","MYBL2","GADD45B","GPX1","PGK1","UBB","VDAC2")
Gon  <- sapply(onco, pick); Gon <- cbind(Gon, ATP5E = pick(c("ATP5E","ATP5F1E")))
gs <- cc[, c("sample","sex","age","tnm_stage","stage_binary","os_time","os_event",
             "kras_raw","braf_raw","tp53_raw","mmr_raw","cimp_raw","CMS","CIT","chemo","risk_score")]
gs <- cbind(gs, as.data.frame(G10), add_cols(G10, Gk2), add_cols(cbind(G10, Gk2), G25),
            add_cols(cbind(G10, Gk2, G25), Gon))
stopifnot(nrow(gs) == 579)
wcsv(gs, "gse39582_analysis.csv")

## ---- 3. xCell scores ----
wcsv(data.frame(cell_type = rownames(xc_t), xc_t, check.names = FALSE), "xcell_tcga_coad.csv")
wcsv(data.frame(cell_type = rownames(xc_g), xc_g, check.names = FALSE), "xcell_gse39582.csv")

## ---- 4. Discovery cohorts: DEG lists and sample sheets ----
wcsv(data.frame(probe = rownames(deg_sig_full), deg_sig_full), "cohort1_DEG_2969.csv")
wcsv(data.frame(probe = rownames(deg_sig), deg_sig), "cohort2_DEG.csv")
wcsv(data.frame(sample = colnames(E1c), dataset = b1, group = grp1), "cohort1_samples.csv")
wcsv(data.frame(sample = colnames(combined_expr), dataset = batch, group = group_all), "cohort2_samples.csv")

## ---- 5. PPI networks (unique undirected edges, gene symbols) ----
edges_sym <- function(df, map) {
  id1 <- pmin(df[[1]], df[[2]]); id2 <- pmax(df[[1]], df[[2]])
  e <- data.frame(string_id1 = id1, string_id2 = id2,
                  gene1 = unname(map[id1]), gene2 = unname(map[id2]), combined_score = df[[3]])
  e <- e[e$string_id1 != e$string_id2, ]
  e[!duplicated(e[, 1:2]), ] }          # unique undirected edges (STRINGdb lists each edge twice)
map1 <- setNames(mp1$gene, mp1$STRING_id); map1 <- map1[!duplicated(names(map1))]
map2 <- setNames(mapped$gene, mapped$STRING_id); map2 <- map2[!duplicated(names(map2))]
wcsv(edges_sym(int1[, c("from","to","combined_score")], map1), "cohort1_ppi_edges.csv")
wcsv(edges_sym(ppi[, c("from","to","combined_score")], map2), "cohort2_ppi_edges.csv")

## ---- 6. Session info ----
writeLines(capture.output(sessionInfo()), file.path(REPO, "sessionInfo.txt"))
cat("\nDone. Files in", out, "\n")
