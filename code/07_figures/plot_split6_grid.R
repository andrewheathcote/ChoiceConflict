## split6 for two models as one A3 landscape page (Andrew, 2026-11-12): 4 x 6, rows LNR_MiSi increasers / decreasers
## then RDM_ViSi increasers / decreasers; columns delta for Simon, Stroop, ANT then CAF for the same tasks, task names
## above the columns. Delta y axis on the left, CAF y axis on the right, each with its outer label as in best5.
## Groups, slopes and panels as plot_split6.R (the groups come from the data, so they are the same people for both
## models); y limits shared per measure across the page from a first automatic pass, CAF capped at 100.
##   Rscript code/07_figures/plot_split6_grid.R [X] [MODEL1 MODEL2]   (default 100 LNR_MiSi RDM_ViSi)
## Writes figures/split6_grid_X<X>.pdf.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
if (LEVEL != "subj") stop("the split needs subject-level predictives: hyper draws are new people (PP_LEVEL=subj)")
source(file.path(CODE, "lib", "predict.R"))                     # stats_one()
a <- commandArgs(trailingOnly = TRUE)
X      <- if (length(a) > 0) as.integer(a[1]) else 100
MODELS <- if (length(a) > 2) a[2:3] else c("LNR_MiSi", "RDM_ViSi")
TASKS  <- c("Simon","Stroop","ANT")
A <- read_tab("adequacy_display.csv")
slope_d <- function(qc, qi) apply(qi - qc, 1, HALF)
slope_c <- function(aI, aC) apply(aI, 1, HALF) - apply(aC, 1, HALF)
grp <- function(s, side) { o <- names(s)[order(s, na.last = NA)]; if (side == "inc") tail(o, X) else head(o, X) }

DAT <- PP <- G <- list()
for (TK in TASKS) {
  for (M in MODELS) { r <- A[A$spec == paste(TK, M, sep = "_"), ][1, ]
    PP[[paste(TK, M)]] <- readRDS(raw_file(r$label)); if (is.null(DAT[[TK]])) DAT[[TK]] <- dat_of(r$label)
    cat(sprintf("%-6s %s %s\n", TK, r$label, r$class)) }
  subs <- levels(factor(DAT[[TK]]$subjects)); ob <- stats_one(DAT[[TK]], subs)
  sl <- list(delta = setNames(slope_d(ob$sqc, ob$sqi), subs), caf = setNames(slope_c(ob$scafI, ob$scafC), subs))
  for (meas in names(sl)) for (side in c("inc", "dec")) G[[paste(TK, meas, side)]] <- grp(sl[[meas]], side)
}
draw <- function(M, TK, meas, side, YL, ...) {
  g <- G[[paste(TK, meas, side)]]; P <- PP[[paste(TK, M)]]
  draw_emc2(DAT[[TK]][DAT[[TK]]$subjects %in% g, ], P[P$subjects %in% g, ], meas, YL, ...)
}
ROWS <- expand.grid(side = c("inc", "dec"), model = MODELS, stringsAsFactors = FALSE)[, 2:1]
tf <- tempfile(fileext = ".pdf"); pdf(tf, width = 16.54, height = 11.69)   # pass 1: automatic limits
par(mfrow = c(4, 6), mar = c(2.2, 2.4, 0.5, 2.4))
LIM <- list()
for (meas in c("delta", "caf")) { v <- c()
  for (i in seq_len(nrow(ROWS))) for (TK in TASKS) v <- range(c(v, draw(ROWS$model[i], TK, meas, ROWS$side[i], NULL)))
  LIM[[meas]] <- v }
invisible(dev.off()); unlink(tf)
LIM$caf[2] <- min(LIM$caf[2], 100)

OUT <- file.path(FIGDIR, emc2_name("split6", sprintf("grid_X%d", X)))
pdf(OUT, width = 16.54, height = 11.69)                                    # A3 landscape
par(mfrow = c(4, 6), oma = c(2.4, 6.4, 1.8, 2.6), mgp = c(1.4, 0.45, 0), tcl = -0.25)
for (i in seq_len(nrow(ROWS))) for (j in 1:6) {
  meas <- if (j <= 3) "delta" else "caf"; TK <- TASKS[(j - 1) %% 3 + 1]
  par(mfg = c(i, j), mar = if (meas == "delta") c(2.2, 2.4, 0.5, 0.4) else c(2.2, 0.4, 0.5, 2.4))
  if (meas == "delta") invisible(draw(ROWS$model[i], TK, meas, ROWS$side[i], LIM[[meas]]))
  else { invisible(draw(ROWS$model[i], TK, meas, ROWS$side[i], LIM[[meas]], yaxt = "n")); axis(4) } }
cols <- (1:6 - 0.5) / 6; rows <- (4:1 - 0.5) / 4                          # column and row centres
mtext(rep(TASKS, 2), side = 3, outer = TRUE, line = 0.3, cex = 1.15, font = 2, at = cols)
## left labels, out from the axis (Andrew, 2026-11-12): Increasers / Decreasers per row, the model per pair of rows,
## then the delta axis label
mtext(c(inc = "Increasers", dec = "Decreasers")[ROWS$side], side = 2, outer = TRUE, line = 0.6, cex = 1, font = 2, at = rows)
mtext(sub("_", " ", MODELS), side = 2, outer = TRUE, line = 2.4, cex = 1.15, font = 2, family = "Palatino",
      at = sapply(MODELS, function(M) mean(rows[ROWS$model == M])))
mtext("RT (I-C) (sec)", side = 2, outer = TRUE, line = 4.4, cex = 0.95)
mtext("CAF (%)", side = 4, outer = TRUE, line = 1.0, cex = 0.95)
mtext(c("Average RT (seconds)", "Bin Centre (%)"), side = 1, outer = TRUE, line = 1.0, cex = 0.9, at = c(0.25, 0.75))
invisible(dev.off())
if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded")   # Palatino is not a base PDF font
cat("wrote", OUT, "\n")
