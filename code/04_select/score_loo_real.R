## Real conflict fits: model rankings by PSIS-LOO vs best-draw DIC (and WAIC), from derived/tables/loo_real.csv
## (code/04_select/loo_real.R) and the pointwise elpd in derived/loo_real/<label>.rds.
##   Rscript code/04_select/score_loo_real.R
## Outputs: derived/tables/loo_real_within.csv  per task x architecture: best model by each criterion, whether LOO and
##            DIC agree, the LOO gap between the two choices with its SE (2 sqrt(n var(pointwise difference)));
##          derived/tables/loo_real_arch.csv    per task: architecture order (best model per architecture) by each
##            criterion, with LOO and DIC differences from the best architecture.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
R <- read_tab("loo_real.csv"); PWD <- file.path(DERIVED, "loo_real")
pw <- function(l) readRDS(file.path(PWD, paste0(l, ".rds")))$elpd_loo
se_d <- function(l1, l2) { if (l1 == l2) return(0); x <- pw(l1); y <- pw(l2); 2 * sqrt(length(x) * var(x - y)) }
## Stroop's diffusion models use accuracy coding (dat_acc), so race and diffusion LOO are compared with that caveat
within <- do.call(rbind, lapply(split(R, list(R$task, R$arch), drop = TRUE), function(x) {
  bL <- which.min(x$LOO); bD <- which.min(x$DIC); bW <- which.min(x$WAIC)
  data.frame(task = x$task[1], arch = x$arch[1], n = nrow(x), best_LOO = x$spec[bL], best_DIC = x$spec[bD], best_WAIC = x$spec[bW],
             agree_LOO_DIC = bL == bD, LOO_of_DIC_choice_minus_best = x$LOO[bD] - x$LOO[bL],
             se = se_d(x$label[bD], x$label[bL]), DIC_of_LOO_choice_minus_best = x$DIC[bL] - x$DIC[bD],
             rank_of_DIC_choice_by_LOO = rank(x$LOO)[bD], k_gt_0.7_max = max(x$k_gt_0.7), row.names = NULL) }))
within <- within[order(within$task, within$arch), ]
write_out(within, tab_path("loo_real_within.csv"))
arch <- do.call(rbind, lapply(split(R, R$task), function(x) {
  best <- function(cr) { b <- do.call(rbind, lapply(split(x, x$arch), function(y) y[which.min(y[[cr]]), ])); b[order(b[[cr]]), ] }
  bl <- best("LOO"); bd <- best("DIC"); bw <- best("WAIC")
  data.frame(task = x$task[1], order_LOO = paste(bl$arch, collapse = " < "), order_DIC = paste(bd$arch, collapse = " < "),
             order_WAIC = paste(bw$arch, collapse = " < "), same_order = identical(bl$arch, bd$arch),
             dLOO = paste(sprintf("%s %.0f", bl$arch, bl$LOO - bl$LOO[1]), collapse = ", "),
             dDIC = paste(sprintf("%s %.0f", bd$arch, bd$DIC - bd$DIC[1]), collapse = ", "), row.names = NULL) }))
write_out(arch, tab_path("loo_real_arch.csv"))
print(within, row.names = FALSE, digits = 4); cat("\n"); print(arch, row.names = FALSE)
