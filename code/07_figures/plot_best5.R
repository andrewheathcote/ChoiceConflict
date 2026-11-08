## best5 drawn with EMC2's plot_delta / plot_caf (see emc2_plot_lib.R): 2 x 5 per task, delta on top, CAF below,
## columns WDM DDM RDM LBA LNR. The same models as plot_best5_display.R: LNR MiSi, RDM ViSi, LBA ViSi fixed; WDM / DDM
## the best for the row's measure (median misfit table, as the paper's display run) among models that met the stop
## rule, are drawn and are not degenerate; Stroop WDM BcVc and DDM BcVc fixed (FIXED_TASK, as plot_best5_display.R).
## y limits fixed and shared across tasks (Andrew, 2026-10-21): delta 0 to .15, CAF 70-100.
##   Rscript code/07_figures/plot_best5.R      (needs derived/pp_raw/<label>_subj.rds from 06_predict/pp_subj.R raw)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
OUTDIR <- FIGDIR
ARCHS  <- c("WDM","DDM","RDM","LBA","LNR")
FIXED  <- c(LNR = "LNR_MiSi", RDM = "RDM_ViSi", LBA = "LBA_ViSi")
## per task, fixed in both rows (Andrew, 2026-10-21): Stroop's diffusion models (accuracy-coded, no z) top out
## at BcVcSc, which has no usable fit, so WDM takes BcVc (best converged by DIC that nests Vc) and DDM Vc
FIXED_TASK <- list(Stroop = c(WDM = "WDM_BcVc", DDM = "DDM_BcVc"))   # DDM Vc -> BcVc once the campaign DDMs finished (Andrew, 2026-11-02)

A  <- read_tab("adequacy_display_median.csv")
DI <- read_tab("display_jm.csv")
MK <- read_tab("mech_key.csv")              # one mechanism key (03_classify/mech_key.R)
A$variant <- sub("^[^_]+_", "", A$spec)
A$mech <- MK$mech[match(A$spec, MK$spec)]                     # mechanisms by spec, never file names
A$drawn <- !(A$arch %in% c("WDM","DDM") & grepl("V", A$mech)) & nchar(A$mech) <= 4
inplay  <- DI$spec[DI$in_play]
pick <- function(tk, ar, meas) {
  if (ar %in% names(FIXED)) { s <- A[A$task == tk & A$variant == FIXED[[ar]], ]; return(s[1, ]) }
  ft <- FIXED_TASK[[tk]]
  if (ar %in% names(ft)) { s <- A[A$task == tk & A$variant == ft[[ar]], ]; return(s[1, ]) }
  s <- A[A$task == tk & A$arch == ar & A$class == "converged" & A$drawn & !is.na(A$rel) & A$degen == 0, ]
  if (!nrow(s)) return(NULL)
  s[which.min(if (meas == "delta") s$rel else s$caf), ]
}

for (TK in c("Simon","Stroop","ANT")) {
  cells <- list()
  for (ar in ARCHS) for (meas in c("delta","caf")) cells[[paste(ar, meas)]] <- list(ar = ar, meas = meas, r = pick(TK, ar, meas))
  if (nzchar(Sys.getenv("LIST_ONLY"))) { for (z in cells) if (!is.null(z$r)) cat("NEED", z$r$label, "\n"); next }   # 07_figures/raw_needed.R
  PP <- list(); DAT <- list()                                # each fit's predictives read once
  for (z in cells) if (!is.null(z$r) && have_raw(z$r$label) && is.null(PP[[z$r$label]])) {
    PP[[z$r$label]] <- readRDS(raw_file(z$r$label)); DAT[[z$r$label]] <- dat_of(z$r$label) }
  draw <- function(z, YL) {
    if (is.null(z$r) || !have_raw(z$r$label)) {
      blank(if (is.null(z$r) && z$ar == "DDM" && any(startsWith(inplay, paste0(TK, "_DDM_")))) "fit running"
            else if (is.null(z$r)) "no model" else "no predictives")
      return(c(NA, NA)) }
    draw_emc2(DAT[[z$r$label]], PP[[z$r$label]], z$meas, YL)
  }
  LIM <- list(delta = c(0, .15), caf = c(70, 100))
  OUT <- file.path(OUTDIR, emc2_name("best5", TK))
  pdf(OUT, width = 15, height = 6.6)
  par(mfcol = c(2, 5), mar = c(2.4, 2.4, 0.6, 0.4), oma = c(2.2, 2.6, 1.6, 0.4), mgp = c(1.4, 0.45, 0), tcl = -0.25)
  for (z in cells) {
    invisible(draw(z, LIM[[z$meas]]))
    if (!is.null(z$r) && have_raw(z$r$label)) label_emc2(paste0(sub("^[A-Z]+_", "", z$r$variant), tag(z$r$class)), z$meas, cex = 0.95 * 1.5)   # 1.5x (Andrew, 4 Nov)
  }
  mtext(ARCHS, side = 3, outer = TRUE, line = 0.25, cex = 1.05, font = 2, at = (seq_along(ARCHS) - 0.5) / length(ARCHS))
  mtext(c("CAF (%)", "RT (I-C) (sec)"), side = 2, outer = TRUE, line = 1.0, cex = 0.9, at = c(0.25, 0.75))
  mtext("Bin Centre (%)  /  Average RT (seconds)", side = 1, outer = TRUE, line = 0.9, cex = 0.85)
  invisible(dev.off()); cat("wrote", OUT, "\n")
  for (z in cells) if (!is.null(z$r)) cat(sprintf("  %-4s %-5s %-24s %s\n", z$ar, z$meas, z$r$label, z$r$class))
  rm(PP, DAT); gc(verbose = FALSE)
}
