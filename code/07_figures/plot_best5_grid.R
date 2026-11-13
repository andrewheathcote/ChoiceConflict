## best5 as one A3 landscape page (Andrew, 2026-11-12): 6 x 5, delta functions in rows 1-3 and CAFs in rows 4-6,
## tasks in the order Simon, Stroop, ANT within each block, columns WDM DDM RDM LBA LNR. Models, panels, labels and
## y limits as plot_best5.R (same pick rule, FIXED / FIXED_TASK, delta 0 to .15, CAF 70-100); only the layout differs.
## Transposed (t, Andrew, 2026-11-12): 5 x 6, rows the model types, columns delta for Simon, Stroop, ANT then CAF for
## the same tasks, task names above the columns; delta y axis on the left, CAF y axis on the right.
## Model labels in Palatino, whose i cannot pass for an l (fonts embedded with Ghostscript).
##   Rscript code/07_figures/plot_best5_grid.R [t]   -> figures/best5_grid.pdf, or best5_grid_5x6.pdf with t
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
OUTDIR <- FIGDIR
ARCHS  <- c("WDM","DDM","RDM","LBA","LNR")
TASKS  <- c("Simon","Stroop","ANT")
FIXED  <- c(LNR = "LNR_MiSi", RDM = "RDM_ViSi", LBA = "LBA_ViSi")
FIXED_TASK <- list(Stroop = c(WDM = "WDM_BcVc", DDM = "DDM_BcVc"))

A  <- read_tab("adequacy_display_median.csv")
DI <- read_tab("display_jm.csv")
MK <- read_tab("mech_key.csv")
A$variant <- sub("^[^_]+_", "", A$spec)
A$mech <- MK$mech[match(A$spec, MK$spec)]
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
LIM <- list(delta = c(0, .15), caf = c(70, 100))
TR <- identical(commandArgs(trailingOnly = TRUE)[1], "t")
FONT <- "Palatino"

OUT <- file.path(OUTDIR, emc2_name("best5", if (TR) "grid_5x6" else "grid"))
pdf(OUT, width = 16.54, height = 11.69)                      # A3 landscape
if (TR) par(mfrow = c(5, 6), oma = c(2.4, 3.6, 1.8, 2.6), mgp = c(1.3, 0.45, 0), tcl = -0.25) else
par(mfrow = c(6, 5), mar = c(2.2, 2.4, 0.5, 0.4), oma = c(2.4, 5.2, 1.8, 0.4), mgp = c(1.3, 0.45, 0), tcl = -0.25)
for (ti in seq_along(TASKS)) {                               # one task's predictives in memory at a time
  TK <- TASKS[ti]; PP <- list(); DAT <- list()
  for (ai in seq_along(ARCHS)) for (meas in c("delta","caf")) {
    ar <- ARCHS[ai]; r <- pick(TK, ar, meas)
    if (!is.null(r) && have_raw(r$label) && is.null(PP[[r$label]])) { PP[[r$label]] <- readRDS(raw_file(r$label)); DAT[[r$label]] <- dat_of(r$label) }
    if (TR) par(mfg = c(ai, if (meas == "delta") ti else 3 + ti),
                mar = if (meas == "delta") c(2.2, 2.4, 0.5, 0.4) else c(2.2, 0.4, 0.5, 2.4))
    else par(mfg = c(if (meas == "delta") ti else 3 + ti, ai))
    if (is.null(r) || !have_raw(r$label)) {
      blank(if (is.null(r) && ar == "DDM" && any(startsWith(inplay, paste0(TK, "_DDM_")))) "fit running"
            else if (is.null(r)) "no model" else "no predictives")
    } else {
      if (TR && meas == "caf") { draw_emc2(DAT[[r$label]], PP[[r$label]], meas, LIM[[meas]], yaxt = "n"); axis(4) }
      else draw_emc2(DAT[[r$label]], PP[[r$label]], meas, LIM[[meas]])
      op <- par(family = if (TR) FONT else ""); label_emc2(paste0(sub("^[A-Z]+_", "", r$variant), tag(r$class)), meas, cex = 1.2); par(op)
    }
    cat(sprintf("  %-6s %-4s %-5s %s\n", TK, ar, meas, if (is.null(r)) "-" else paste(r$label, r$class)))
  }
  rm(PP, DAT); gc(verbose = FALSE)
}
if (TR) {
  mtext(rep(TASKS, 2), side = 3, outer = TRUE, line = 0.3, cex = 1.15, font = 2, at = (1:6 - 0.5) / 6)
  mtext(ARCHS, side = 2, outer = TRUE, line = 2.2, cex = 1.15, font = 2, at = (5:1 - 0.5) / 5)
  mtext("RT (I-C) (sec)", side = 2, outer = TRUE, line = 0.6, cex = 0.95)
  mtext("CAF (%)", side = 4, outer = TRUE, line = 1.0, cex = 0.95)
  mtext(c("Average RT (seconds)", "Bin Centre (%)"), side = 1, outer = TRUE, line = 1.0, cex = 0.9, at = c(0.25, 0.75))
} else {
mtext(ARCHS, side = 3, outer = TRUE, line = 0.3, cex = 1.15, font = 2, at = (seq_along(ARCHS) - 0.5) / length(ARCHS))
rows <- (6:1 - 0.5) / 6                                      # row centres, top to bottom
mtext(rep(TASKS, 2), side = 2, outer = TRUE, line = 3.4, cex = 1.05, font = 2, at = rows)
mtext(c("RT (I-C) (sec)", "CAF (%)"), side = 2, outer = TRUE, line = 1.4, cex = 0.95, at = c(0.75, 0.25))
mtext("Average RT (seconds) (delta, rows 1-3)  /  Bin Centre (%) (CAF, rows 4-6)", side = 1, outer = TRUE, line = 1.0, cex = 0.9)
}
invisible(dev.off())
if (TR) { if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded") }
cat("wrote", OUT, "\n")
