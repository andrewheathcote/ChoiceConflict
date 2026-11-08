## models6 drawn with EMC2's plot_delta / plot_caf (see emc2_plot_lib.R): per task 2 x 6, delta on top, CAF below,
## columns RDM VciSci / ViSi, RDM BciSci / BiSi, LNR MciSci / MiSi, each from its displayed fit (as
## plot_models6_display.R). y limits fixed and shared across tasks (Andrew, 2026-10-21): delta -.0025 to .085,
## CAF 70-100.
##   Rscript code/07_figures/plot_models6.R    (needs derived/pp_raw/<label>_subj.rds from 06_predict/pp_subj.R raw)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
OUTDIR <- FIGDIR
MODELS <- c("RDM_VciSci","RDM_ViSi","RDM_BciSci","RDM_BiSi","LNR_MciSci","LNR_MiSi")
A <- read_tab("adequacy_display.csv")
YLIM <- list(delta = c(-.0025, .085), caf = c(70, 100))
LABEL_BOTTOM <- character(0)          # Simon_RDM_BciSci's delta label back at the top (Andrew, 4 Nov; moved down 2026-10-21)
LABEL_CEX <- 0.95 * 1.5               # labels 1.5x (Andrew, 4 Nov)
for (TK in c("Simon","Stroop","ANT")) {
  OUT <- file.path(OUTDIR, emc2_name("models6", TK))
  pdf(OUT, width = 16, height = 6.4)
  par(mfcol = c(2, 6), mar = c(2.4, 2.4, 0.6, 0.4), oma = c(2.2, 2.8, 0.4, 0.4), mgp = c(1.4, 0.45, 0), tcl = -0.25)
  for (v in MODELS) {
    r <- A[A$spec == paste(TK, v, sep = "_"), ][1, ]
    if (!have_raw(r$label)) { blank("no predictives"); blank("no predictives"); next }
    dat <- dat_of(r$label); pp <- readRDS(raw_file(r$label))
    for (meas in c("delta", "caf")) {
      invisible(draw_emc2(dat, pp, meas, YLIM[[meas]]))
      label_emc2(sub("^[A-Z]+_", "", v), meas,
                 bottom = r$spec %in% LABEL_BOTTOM, cex = LABEL_CEX)   # variant only: no convergence tags (Andrew, 4 Nov)
    }
    rm(pp, dat); gc(verbose = FALSE)
  }
  mtext(c("CAF (%)", "RT (I-C) (sec)"), side = 2, outer = TRUE, line = 1.2, cex = 0.9, at = c(0.25, 0.75))
  mtext("Bin Centre (%)  /  Average RT (seconds)", side = 1, outer = TRUE, line = 0.9, cex = 0.85)
  invisible(dev.off()); cat("wrote", OUT, "\n")
}
