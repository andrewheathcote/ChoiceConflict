## Tables of the supplement, generated from the archive's tables (so the numbers cannot drift from them).
##   Rscript doc/make_tables.R      -> doc/src/tables/*.md (Markdown pipe tables, included by doc/build.sh)
CODE <- file.path(dirname(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))), "code")
source(file.path(CODE, "config.R"))
OUT <- file.path(ROOT, "doc", "src", "tables"); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
md <- function(d, file, caption) {
  fmt <- function(x) if (is.numeric(x)) vapply(x, function(v) if (is.na(v)) "" else formatC(v, format = "f", digits = if (v != round(v)) 3 else 0, big.mark = ","), "") else ifelse(is.na(x), "", as.character(x))
  d[] <- lapply(d, fmt)
  l <- c(caption, "", paste0("| ", paste(names(d), collapse = " | "), " |"), paste0("|", paste(rep("---", ncol(d)), collapse = "|"), "|"),
         apply(d, 1, function(r) paste0("| ", paste(r, collapse = " | "), " |")))
  writeLines(l, file.path(OUT, file))
}
## S3.1 direction of the unrestricted congruency terms (before the Q6 refits for Simon LNR)
q <- read_tab("q6_direction_before.csv")
q <- data.frame(fit = q$label, term = q$par, `group mean (median)` = round(q$mu_median, 3), `P(group mean > 0)` = round(q$P_mu_pos, 3),
                `subjects with posterior median > 0` = round(q$subj_median_pos, 3), `subjects P(> 0) > .95` = round(q$subj_P95_pos, 3),
                `subjects P(< 0) > .95` = round(q$subj_P95_neg, 3), check.names = FALSE)
md(q, "S3_1.md", "**Table S3.1.** Posterior direction of the congruency terms that were estimated without an order restriction: Simon LNR m_CIincong (the original fits, before the Q6 refits restricted it) and Stroop WDM/DDM v_CIcong. The predicted direction is positive.")
## S6.1 classes by task and architecture
D <- read_tab("display_jm.csv")
lev <- c("converged", "ridge-only, shown as converged", "converged, subject-mode caveat", "ridge-only, Rhat >= 1.5", "not converged")
t6 <- as.data.frame.matrix(table(paste(D$task, D$arch), factor(D$class, lev)))
t6 <- cbind(`task / architecture` = rownames(t6), t6, models = rowSums(t6), displayed = as.vector(tapply(D$display_converged, paste(D$task, D$arch), sum)[rownames(t6)]))
t6 <- t6[order(match(sub(" .*", "", t6[[1]]), TASKS), sub(".* ", "", t6[[1]])), ]
t6 <- rbind(t6, cbind(`task / architecture` = "total", as.data.frame(t(colSums(t6[, -1])))))
md(t6, "S6_1.md", "**Table S6.1.** Convergence class of each model's display run (S6). Displayed as converged: converged plus ridge-only shown as converged.")
## S7.1 simple choice: DIC, WAIC, PSIS-LOO and weights for each family's best fit (04_select/compare_choice.R)
cc <- read_tab("compare_choice.csv")
d7 <- function(x) round(x - min(x))
cc <- data.frame(family = cc$model, fit = cc$label, `DIC difference` = d7(cc$DIC), `WAIC difference` = d7(cc$WAIC),
                 `LOO difference` = d7(cc$LOO), `pseudo-BMA+ weight` = round(cc$wLOO, 3), `stacking weight` = round(cc$wStack, 3),
                 `observations with k > 0.7` = cc$k_bad, check.names = FALSE)
md(cc[order(cc$`LOO difference`), ], "S7_1.md", "**Table S7.1.** Simple choice: the best fit of each family compared by best-draw conditional DIC, WAIC and PSIS-LOO (differences from the best on the deviance scale), with pseudo-BMA+ (Bayesian bootstrap) and stacking weights from PSIS-LOO, and the number of observations (of 77,133) with Pareto k > 0.7.")
## S9.1 recovery study 5 per cell
c5 <- read.csv(file.path(ROOT, "recovery", "conflict", "recovery_rec5_cells.csv"), stringsAsFactors = FALSE)
c5 <- data.frame(cell = paste(c5$task, c5$arch), `G1 (one)` = c5$G1, `G2 (two)` = c5$G2, `G1 recovered` = sprintf("%d/%d", c5$recovered_g1, c5$sets_g1),
                 `false alarms (max margin)` = sprintf("%d (%.0f)", c5$false_alarms, c5$fa_max_margin),
                 `G2 recovered` = sprintf("%d/%d", c5$recovered_g2, c5$sets_g2), `misses (max margin)` = sprintf("%d (%.0f)", c5$misses, c5$miss_max_margin),
                 `real data: best two-mechanism - G1` = round(c5$real_best2_minus_g1), check.names = FALSE)
md(c5, "S9_1.md", "**Table S9.1.** Recovery study 5 by cell: 5 data sets per generator, 4 candidates each. False alarm: a two-mechanism model wins on G1 data; miss: G1 wins on G2 data (margins in DIC). Last column: the real data's DIC difference (negative favours two mechanisms).")
cat("tables written to", OUT, "\n")
