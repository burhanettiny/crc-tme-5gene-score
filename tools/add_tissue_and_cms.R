## ============================================================
## tools/add_tissue_and_cms.R   (revision fix, run once)
## GSE39582 contains 566 primary tumours + 19 non-tumoural mucosa samples
## (Marisa et al. 2013). The first export kept 17 mucosa samples that carry
## patient-level survival data. This script
##   (1) adds a `tissue` column to data/processed/gse39582_analysis.csv, and
##   (2) recomputes CMS labels (CMScaller) on primary tumours only.
## Run from the repository root.
## ============================================================
suppressPackageStartupMessages({ library(GEOquery); library(Biobase); library(CMScaller) })
SERIES <- "D:/1-Kisisel/1-Akademik/5-Yayinlar/29-Col-Can/GSE39582_series_matrix.txt.gz"

g  <- getGEO(filename = SERIES)              # also fetches the GPL570 annotation
p  <- pData(g); fd <- fData(g)
gs <- read.csv("data/processed/gse39582_analysis.csv", check.names = FALSE)

## ---- (1) tissue ----
gs$tissue <- ifelse(grepl("Adenocarcinoma", p[gs$sample, "source_name_ch1"]), "tumor", "normal")
print(table(gs$tissue))

## ---- (2) CMS on primary tumours only ----
tum  <- gs$sample[gs$tissue == "tumor"]
ecol <- intersect(c("ENTREZ_GENE_ID", "Gene ID"), names(fd))[1]
ent  <- as.character(fd[[ecol]])
E    <- exprs(g)[, tum]
ok   <- !is.na(ent) & ent != "" & !grepl("///", ent)
E    <- E[ok, ]; ent <- ent[ok]
vv   <- apply(E, 1, var)
keep <- tapply(seq_along(ent), ent, function(i) i[which.max(vv[i])])
E    <- E[unlist(keep), ]; rownames(E) <- names(keep)
set.seed(123)
res  <- CMScaller(t(scale(t(E))), rowNames = "entrez", RNAseq = FALSE, nPerm = 1000, FDR = 0.05, doPlot = FALSE)
lab  <- as.character(res$prediction); lab[is.na(lab)] <- "Unclassified"
gs$CMS <- NA
gs$CMS[match(tum, gs$sample)] <- lab[match(tum, rownames(res))]
print(table(gs$CMS[gs$tissue == "tumor"], useNA = "ifany"))

write.csv(gs, "data/processed/gse39582_analysis.csv", row.names = FALSE)
cat("Updated data/processed/gse39582_analysis.csv (tumours with survival data:",
    sum(gs$tissue == "tumor"), ")\n")
