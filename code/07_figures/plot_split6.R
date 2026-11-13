## split6 drawn with EMC2's plot_delta / plot_caf (see emc2_plot_lib.R), subject level only, as the archive's
## plot_figB.R: 2 x 6 for one model (RDM_ViSi), rows increasers / decreasers (top / bottom X real subjects), columns
## delta for the three tasks then CAF for the three tasks; each column splits on the measure it shows. The slopes are
## those of plot_split6_display.R (data, mean of points 4-5 minus mean of points 1-2, rank_subjects.R's rule), so the
## groups are the same people. y limits shared per measure from a first automatic pass, CAF capped at 100.
##   Rscript code/07_figures/plot_split6.R [X] [MODEL]
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
if (LEVEL != "subj") stop("the split needs subject-level predictives: hyper draws are new people (PP_LEVEL=subj)")
source(file.path(CODE, "lib", "predict.R"))                     # stats_one()
a <- commandArgs(trailingOnly = TRUE)
X     <- if (length(a) > 0) as.integer(a[1]) else 100
MODEL <- if (length(a) > 1) a[2] else "RDM_ViSi"
TASKS <- c("Simon","Stroop","ANT")
OUTDIR <- FIGDIR
A <- read_tab("adequacy_display.csv")
slope_d <- function(qc, qi) apply(qi - qc, 1, HALF)
slope_c <- function(aI, aC) apply(aI, 1, HALF) - apply(aC, 1, HALF)
grp <- function(s, side) { o <- names(s)[order(s, na.last = NA)]; if (side == "inc") tail(o, X) else head(o, X) }

DAT <- PP <- G <- list()
for (TK in TASKS) {
  r <- A[A$spec == paste(TK, MODEL, sep = "_"), ][1, ]
  DAT[[TK]] <- dat_of(r$label); PP[[TK]] <- readRDS(raw_file(r$label))
  subs <- levels(factor(DAT[[TK]]$subjects)); ob <- stats_one(DAT[[TK]], subs)
  sl <- list(delta = setNames(slope_d(ob$sqc, ob$sqi), subs), caf = setNames(slope_c(ob$scafI, ob$scafC), subs))
  for (meas in names(sl)) for (side in c("inc", "dec")) G[[paste(TK, meas, side)]] <- grp(sl[[meas]], side)
  cat(sprintf("%-6s %s %s\n", TK, r$label, r$class))
}
draw <- function(TK, meas, side, YL) {
  g <- G[[paste(TK, meas, side)]]
  draw_emc2(DAT[[TK]][DAT[[TK]]$subjects %in% g, ], PP[[TK]][PP[[TK]]$subjects %in% g, ], meas, YL)
}
tf <- tempfile(fileext = ".pdf"); pdf(tf, width = 16, height = 6.4)       # pass 1: automatic limits
par(mfrow = c(2, 6), mar = c(2.4, 2.4, 0.6, 0.4), oma = c(2.2, 2.8, 0.4, 0.4))
LIM <- list()
for (meas in c("delta", "caf")) { v <- c()
  for (side in c("inc", "dec")) for (TK in TASKS) v <- range(c(v, draw(TK, meas, side, NULL)))
  LIM[[meas]] <- v }
invisible(dev.off()); unlink(tf)
LIM$caf[2] <- min(LIM$caf[2], 100)
OUT <- file.path(OUTDIR, emc2_name("split6", sprintf("%s_X%d", MODEL, X)))
pdf(OUT, width = 16, height = 6.4)
par(mfrow = c(2, 6), mar = c(2.4, 2.4, 0.6, 0.4), oma = c(2.2, 2.8, 0.4, 0.4), mgp = c(1.4, 0.45, 0), tcl = -0.25)
for (side in c("inc", "dec")) for (meas in c("delta", "caf")) for (TK in TASKS) {
  invisible(draw(TK, meas, side, LIM[[meas]])); label_emc2(TK, meas, cex = 1.5 * 0.95) }   # task labels 1.5x (Andrew, 2026-11-10)
mtext(c("Decreasers", "Increasers"), side = 2, outer = TRUE, line = 1.2, cex = 0.95, font = 2, at = c(0.25, 0.75))
mtext("Average RT (seconds)   /   Bin Centre (%)", side = 1, outer = TRUE, line = 0.9, cex = 0.85)
invisible(dev.off()); cat("wrote", OUT, "\n")
