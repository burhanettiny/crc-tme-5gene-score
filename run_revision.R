## Runs all revision analyses (scripts 08-13) from data/processed/.
## Each script is run in a fresh environment; outputs go to results/revision/.
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)

scripts <- c("R/08_revision_model_diagnostics.R", "R/09_revision_platform.R",
             "R/10_revision_subgroups_utility.R", "R/11_revision_signatures.R",
             "R/12_revision_xcell.R", "R/13_revision_ppi_batch.R")
for (s in scripts) {
  cat("\n\n#################", s, "#################\n")
  sink(file.path("results/revision", sub("\\.R$", ".log", basename(s))), split = TRUE)
  tryCatch(source(s, local = new.env()), error = function(e) message("ERROR in ", s, ": ", conditionMessage(e)))
  sink()
}
writeLines(capture.output(sessionInfo()), "results/revision/sessionInfo.txt")
