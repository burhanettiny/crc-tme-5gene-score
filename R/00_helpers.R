## ============================================================
## 00_helpers.R — shared data loading and helper functions
## Sourced by scripts 08–13. Run all scripts from the repository root.
## ============================================================
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(timeROC)
})
dir.create("results/revision", recursive = TRUE, showWarnings = FALSE)
outp <- function(f) file.path("results/revision", f)

TME10 <- c("IL6","IL1B","CXCL8","MYC","PTGS2","COL1A1","IL1A","TIMP1","CXCL12","MMP3")
G5    <- c("IL6","CXCL8","PTGS2","TIMP1","MMP3")
K2    <- c("MYC","CDK1","CCNB1","TOP2A","CD44","CCNA2","AURKA","CHEK1","MAD2L1","KIF23")

## ---- Processed analysis tables (see tools/export_processed_data.R) ----
tc <- read.csv("data/processed/tcga_coad_analysis.csv", check.names = FALSE)
gs <- read.csv("data/processed/gse39582_analysis.csv",  check.names = FALSE)
cf <- read.csv("data/processed/tcga_lasso_coefficients.csv")
b_orig <- setNames(cf$coef, cf$gene)                       # original LASSO coefficients (10 genes)

tc$stage_f <- factor(tc$stage_binary, levels = c("I-II","III-IV"))
gs$stage_f <- factor(gs$stage_binary, levels = c("I-II","III-IV"))
tc$risk_score <- as.numeric(as.matrix(tc[, names(b_orig)]) %*% b_orig)
stopifnot(isTRUE(all.equal(gs$risk_score,
                           as.numeric(as.matrix(gs[, names(b_orig)]) %*% b_orig))))
y_tcga <- Surv(tc$OS.time, tc$OS)                          # days
X10 <- as.matrix(tc[, TME10]); G10 <- as.matrix(gs[, TME10])

## ---- Transformations ----
zs  <- function(M) scale(M)
int <- function(M) apply(M, 2, function(x) qnorm((rank(x) - 0.5) / length(x)))   # rank-based INT

## ---- LASSO-Cox with fixed seed; Breslow ties locked for glmnet >= 5.0 ----
cvg <- function(X, y, ...) {
  if (packageVersion("glmnet") >= "5.0")
    cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10, cox.ties = "breslow", ...)
  else cv.glmnet(X, y, family = "cox", alpha = 1, nfolds = 10, ...)
}
lasso_fit <- function(X, y, seed = 123) {
  set.seed(seed); cv <- suppressWarnings(cvg(X, y))
  list(cv = cv, b = setNames(as.numeric(coef(cv, s = "lambda.min")), colnames(X)))
}

## ---- Standard evaluation of a continuous score ----
## HR per SD (univariable and stage-adjusted), time-dependent AUC at t_auc[2],
## tertile log-rank test. t_auc in the time unit of `time`.
evaluate_score <- function(score, time, event, stage, t_auc, label) {
  tryCatch({
    d <- data.frame(s = as.numeric(score), time, event,
                    stage = factor(stage, levels = c("I-II","III-IV")))
    d$s_sd <- as.numeric(scale(d$s))
    su <- summary(coxph(Surv(time, event) ~ s_sd, data = d))
    sm <- summary(coxph(Surv(time, event) ~ s_sd + stage, data = d))
    r  <- timeROC(T = d$time, delta = d$event, marker = d$s, cause = 1, times = t_auc, iid = TRUE)
    ci <- confint(r)$CI_AUC
    tq <- quantile(d$s, c(1/3, 2/3)); d$tt <- cut(d$s, c(-Inf, tq, Inf))
    lr <- survdiff(Surv(time, event) ~ tt, data = d)
    data.frame(Score = label, HR_SD = su$conf.int[1,1], L95 = su$conf.int[1,3], U95 = su$conf.int[1,4],
               p = su$coefficients[1,5], HR_SD_stage = sm$conf.int[1,1], p_stage = sm$coefficients[1,5],
               AUC3 = unname(r$AUC[2]), AUC_L = ci[2,1] / 100, AUC_U = ci[2,2] / 100,
               p_logrank = 1 - pchisq(lr$chisq, 2))
  }, error = function(e) { message(label, " -> ", conditionMessage(e)); NULL })
}
TT_TCGA <- c(365, 1095, 1825)   # days
TT_GSE  <- c(12, 36, 60)        # months
