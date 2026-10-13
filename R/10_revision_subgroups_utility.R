## ============================================================
## 10_revision_subgroups_utility.R
## Reviewer comments 3a–b, 4, minor 4, statistical review:
## interaction tests, subgroup table (Table 3), 5-year survival with
## 95% CI, incremental value over stage (C-index, time-dependent AUC,
## DeLong, NRI/IDI, decision curve), adjuvant chemotherapy adjustment.
## Outputs: Table3_subgroups.csv, TableS9–S11, TableS16, FigS6, FigS7
## ============================================================
source("R/00_helpers.R")
suppressPackageStartupMessages({ library(ggplot2); library(survIDINRI); library(dcurves); library(pROC) })
cc <- gs

## ---- 1. Interaction tests (likelihood-ratio) ----
mol <- list(BRAF = c("braf_raw","WT","M"), KRAS = c("kras_raw","WT","M"), TP53 = c("tp53_raw","WT","M"),
            MMR  = c("mmr_raw","pMMR","dMMR"), CIMP = c("cimp_raw","-","+"))
int_tab <- do.call(rbind, lapply(names(mol), function(k) {
  m <- mol[[k]]; d <- cc[!is.na(cc[[m[1]]]) & cc[[m[1]]] %in% m[2:3], ]
  d$grp <- factor(d[[m[1]]], levels = m[2:3])
  m0  <- coxph(Surv(os_time, os_event) ~ risk_score + grp, data = d)
  m1  <- coxph(Surv(os_time, os_event) ~ risk_score * grp, data = d)
  m0s <- coxph(Surv(os_time, os_event) ~ risk_score + grp + stage_f, data = d)
  m1s <- coxph(Surv(os_time, os_event) ~ risk_score * grp + stage_f, data = d)
  ci  <- summary(m1)$conf.int[3, ]
  data.frame(marker = k, n = nrow(d), events = sum(d$os_event), HR_interaction = ci[1], L95 = ci[3], U95 = ci[4],
             p_int = anova(m0, m1)$`Pr(>|Chi|)`[2], p_int_stage = anova(m0s, m1s)$`Pr(>|Chi|)`[2]) }))
int_tab$p_int_FDR <- p.adjust(int_tab$p_int, "BH")
print(int_tab, digits = 3, row.names = FALSE)
write.csv(int_tab, outp("TableS9_interactions.csv"), row.names = FALSE)

## ---- 2. Subgroup table (per-unit HR; FDR across all 10 subgroups jointly) ----
sub_row <- function(d, grp, marker) {
  s  <- summary(coxph(Surv(os_time, os_event) ~ risk_score, data = d))
  sa <- summary(coxph(Surv(os_time, os_event) ~ risk_score + stage_f, data = d))
  tq <- quantile(d$risk_score, c(1/3, 2/3)); d$tt <- cut(d$risk_score, c(-Inf, tq, Inf))
  lr <- survdiff(Surv(os_time, os_event) ~ tt, data = d)
  data.frame(marker = marker, subgroup = grp, n = nrow(d), events = sum(d$os_event),
             HR = s$conf.int[1,1], L95 = s$conf.int[1,3], U95 = s$conf.int[1,4], p = s$coefficients[1,5],
             HR_stage = sa$conf.int[1,1], p_stage = sa$coefficients[1,5], p_logrank = 1 - pchisq(lr$chisq, 2)) }
sub_tab <- sub_row(cc, "Overall", "Overall")
for (k in names(mol)) for (lv in mol[[k]][2:3]) sub_tab <- rbind(sub_tab, sub_row(cc[which(cc[[mol[[k]][1]]] == lv), ], lv, k))
sub_tab$p_FDR <- NA; sub_tab$p_FDR[-1] <- p.adjust(sub_tab$p[-1], "BH")
sub_tab$p_int <- int_tab$p_int[match(sub_tab$marker, int_tab$marker)]
print(sub_tab, digits = 3, row.names = FALSE)
write.csv(sub_tab, outp("Table3_subgroups.csv"), row.names = FALSE)

fp <- sub_tab
fp$label <- factor(sprintf("%s (n=%d)", ifelse(fp$marker == "Overall", "Overall", paste(fp$marker, fp$subgroup)), fp$n),
                   levels = rev(sprintf("%s (n=%d)", ifelse(fp$marker == "Overall", "Overall", paste(fp$marker, fp$subgroup)), fp$n)))
fp$pint <- ""; for (k in names(mol)) { i <- which(fp$marker == k)[1]; fp$pint[i] <- sprintf("P-int = %.2f", fp$p_int[i]) }
g <- ggplot(fp, aes(HR, label)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = L95, xmax = U95)) + geom_point(aes(size = n), shape = 15) +
  geom_text(aes(x = max(U95) * 1.6, label = pint), hjust = 0, size = 3) +
  scale_x_log10(expand = expansion(mult = c(0.05, 0.6))) + scale_size(range = c(1.5, 4), guide = "none") +
  labs(x = "Hazard ratio per unit risk score (95% CI, log scale)", y = NULL) + theme_bw(base_size = 10)
ggsave(outp("FigS6_subgroup_forest.png"), g, width = 7, height = 4.5, dpi = 300)

## ---- 3. 5-year overall survival by tertile (95% CI) + trend test ----
s5 <- function(d, label) {
  tq <- quantile(d$risk_score, c(1/3, 2/3))
  d$tt <- cut(d$risk_score, c(-Inf, tq, Inf), labels = c("Low","Intermediate","High"))
  sm <- summary(survfit(Surv(os_time, os_event) ~ tt, data = d), times = 60)
  pt <- summary(coxph(Surv(os_time, os_event) ~ as.numeric(tt), data = d))$coefficients[1, 5]
  data.frame(cohort = label, tertile = sub("tt=", "", sm$strata), n_at_risk = sm$n.risk,
             OS5 = round(100 * sm$surv, 1), L95 = round(100 * sm$lower, 1), U95 = round(100 * sm$upper, 1),
             p_trend = signif(pt, 3)) }
s5_tab <- rbind(s5(cc, "GSE39582 overall"), s5(cc[which(cc$braf_raw == "WT"), ], "BRAF wild-type"))
print(s5_tab, row.names = FALSE)
write.csv(s5_tab, outp("TableS10_5yr_survival.csv"), row.names = FALSE)

## ---- 4. Incremental value over AJCC stage (GSE39582) ----
d4 <- cc[!is.na(cc$stage_f), ]
m0 <- coxph(Surv(os_time, os_event) ~ stage_f, data = d4, x = TRUE)
m1 <- coxph(Surv(os_time, os_event) ~ stage_f + risk_score, data = d4, x = TRUE)
d4$lp0 <- predict(m0, type = "lp"); d4$lp1 <- predict(m1, type = "lp")
cb <- concordance(m0, m1); dC <- diff(cb$concordance)
seC <- sqrt(sum(diag(cb$var)) - 2 * cb$var[1, 2])
r0 <- timeROC(d4$os_time, d4$os_event, d4$lp0, cause = 1, times = TT_GSE, iid = TRUE)
r1 <- timeROC(d4$os_time, d4$os_event, d4$lp1, cause = 1, times = TT_GSE, iid = TRUE)
p_auc <- timeROC::compare(r0, r1)$p_values_AUC
b36 <- subset(d4, os_event == 1 | os_time >= 36); b36$y36 <- as.integer(b36$os_event == 1 & b36$os_time <= 36)
roc0 <- roc(b36$y36, b36$lp0, quiet = TRUE); roc1 <- roc(b36$y36, b36$lp1, quiet = TRUE)
dl <- roc.test(roc0, roc1, method = "delong")
set.seed(123)
idi <- IDI.INF(indata = d4[, c("os_time","os_event")],
               covs0 = model.matrix(~ stage_f, d4)[, -1, drop = FALSE],
               covs1 = model.matrix(~ stage_f + risk_score, d4)[, -1], t0 = 36, npert = 500)
IDI.INF.OUT(idi)
idi_out <- rbind(idi$m1, idi$m2, idi$m3)          # rows: IDI, continuous NRI, median improvement
inc <- data.frame(
  measure = c("LRT p", "C-index stage", "C-index stage+score", "Delta C (95% CI)", "Delta C p",
              paste0("AUC ", TT_GSE, "m stage"), paste0("AUC ", TT_GSE, "m stage+score"), paste0("AUC ", TT_GSE, "m p"),
              "DeLong 36m AUC stage", "DeLong 36m AUC stage+score", "DeLong p", "DeLong n / events",
              "IDI (95% CI) p", "Continuous NRI (95% CI) p", "Median improvement (95% CI) p"),
  value = c(signif(anova(m0, m1)$`Pr(>|Chi|)`[2], 3), round(cb$concordance[1], 3), round(cb$concordance[2], 3),
            sprintf("%.3f (%.3f to %.3f)", dC, dC - 1.96 * seC, dC + 1.96 * seC), signif(2 * pnorm(-abs(dC / seC)), 3),
            round(r0$AUC, 3), round(r1$AUC, 3), signif(p_auc, 3),
            round(auc(roc0), 3), round(auc(roc1), 3), signif(dl$p.value, 3), paste(nrow(b36), "/", sum(b36$y36)),
            apply(idi_out, 1, function(r) sprintf("%.3f (%.3f to %.3f) p=%.3f", r[1], r[2], r[3], r[4]))))
print(inc, row.names = FALSE)
write.csv(inc, outp("TableS11_incremental_value.csv"), row.names = FALSE)

risk36 <- function(m) 1 - as.numeric(summary(survfit(m, newdata = d4), times = 36)$surv)
d4$Stage_only <- risk36(m0); d4$Stage_plus_5gene <- risk36(m1)
dc <- dca(Surv(os_time, os_event) ~ Stage_only + Stage_plus_5gene, data = d4, time = 36,
          thresholds = seq(0.05, 0.50, by = 0.01))
ggsave(outp("FigS7_DCA_36m.png"), plot(dc, smooth = TRUE) + theme_bw(base_size = 10), width = 6.5, height = 4.5, dpi = 300)
write.csv(as.data.frame(dc$dca)[, c("label","threshold","net_benefit")], outp("TableS11_DCA.csv"), row.names = FALSE)

## ---- 5. Adjustment for adjuvant chemotherapy ----
dk <- cc[!is.na(cc$chemo) & cc$chemo %in% c("Y","N"), ]; dk$chemo <- factor(dk$chemo, levels = c("N","Y"))
mk <- summary(coxph(Surv(os_time, os_event) ~ risk_score + stage_f + chemo, data = dk))
d23 <- dk[dk$tnm_stage %in% c(2, 3, "2", "3"), ]
mk23 <- summary(coxph(Surv(os_time, os_event) ~ risk_score + stage_f + chemo, data = d23))
chemo_tab <- data.frame(set = c(rep("all", 3), rep("stage II-III", 3)), term = c(rownames(mk$conf.int), rownames(mk23$conf.int)),
                        n = c(rep(nrow(dk), 3), rep(nrow(d23), 3)),
                        HR = round(c(mk$conf.int[, 1], mk23$conf.int[, 1]), 3), L95 = round(c(mk$conf.int[, 3], mk23$conf.int[, 3]), 3),
                        U95 = round(c(mk$conf.int[, 4], mk23$conf.int[, 4]), 3), p = signif(c(mk$coefficients[, 5], mk23$coefficients[, 5]), 3))
print(chemo_tab, row.names = FALSE)
write.csv(chemo_tab, outp("TableS16_adjuvant_chemotherapy.csv"), row.names = FALSE)
