## Misfit statistics for every model in display_jm.csv, from its SUBJECT-LEVEL posterior predictive
## (derived/pp_subj/<label>.rds, 06_predict/pp_subj.R: 100 draws, each real subject simulated from its own alpha draw).
## The group curve aggregates the subjects' quantiles: PP_AGG=mean (Vincentised, default) or median; accuracy is
## always the mean over subjects (lib/pp_display_lib.R).
##   delta_ms  RMS over the five quantiles of (median predicted delta - observed delta), ms
##   rel       100 * delta_ms / effect size (Simon 42, Stroop 122, ANT 70 ms)
##   cafI_pp   RMS over the five RT bins of (median predicted - observed accuracy), incongruent, % points
##   cafC_pp   same, congruent;   caf = (cafI_pp + cafC_pp) / 2, what the heat grids plot
##   degen     % of draws whose largest |delta| exceeds 10 s (> 0 excludes the model)
## From newsampler/figs/build_adequacy_display.R, which read the hyper-level primitives in the fit summaries.
##   [PP_AGG=median] Rscript code/06_predict/build_adequacy.R   -> derived/tables/adequacy_display[_median].csv
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "pp_display_lib.R"))
EFFECT <- c(Simon = 42, Stroop = 122, ANT = 70)
D <- read_tab("display_jm.csv")
stat1 <- function(z, tk) {
  pp <- z$pp; obs <- z$obs
  if (is.null(pp) || inherits(pp, "try-error") || is.null(obs)) return(NULL)
  P <- grp_pp(pp); O <- grp_obs(obs)
  d <- P$qi - P$qc
  dm <- sqrt(mean((1000 * (apply(d, 2, median, na.rm = TRUE) - (O$qi - O$qc)))^2))
  data.frame(degen    = 100 * mean(apply(abs(d), 1, max) > 10, na.rm = TRUE),
             delta_ms = dm, rel = 100 * dm / EFFECT[[tk]],
             cafI_pp  = 100 * sqrt(mean((apply(P$aI, 2, median, na.rm = TRUE) - O$aI)^2)),
             cafC_pp  = 100 * sqrt(mean((apply(P$aC, 2, median, na.rm = TRUE) - O$aC)^2)),
             n_bad = pp$n_bad, n_retry = pp$n_retry)
}
rows <- lapply(seq_len(nrow(D)), function(i) {
  dest <- paste0(D$label[i], ".rds")
  s <- if (!is.na(D$label[i]) && have_pp(dest)) stat1(read_summ(dest), D$task[i]) else NULL
  if (is.null(s)) s <- data.frame(degen = NA, delta_ms = NA, rel = NA, cafI_pp = NA, cafC_pp = NA, n_bad = NA, n_retry = NA)
  cbind(D[i, c("spec", "task", "arch", "label", "class", "display_converged", "DIC")], s)
})
A <- do.call(rbind, rows); A$caf <- (A$cafI_pp + A$cafC_pp) / 2; A$level <- "subj"
A <- A[order(A$task, A$arch, A$spec), ]
write_out(A, tab_path(sprintf("adequacy_display%s.csv", aggsuf)))
cat("models:", nrow(A), " with statistics:", sum(!is.na(A$rel)), " displayed without statistics:",
    sum(A$display_converged & is.na(A$rel)), " degenerate:", sum(A$degen > 0, na.rm = TRUE), "\n")
print(A[A$display_converged & is.na(A$rel), c("spec", "label", "class")], row.names = FALSE)
print(A[!is.na(A$degen) & A$degen > 0, c("spec", "label", "class", "degen")], row.names = FALSE)
