## models6 across two tasks (Andrew, 2026-11-12): 2 x 6, the four-mechanism models in row 1 and their two-mechanism
## counterparts in row 2 (LNR MciSci / MiSi, RDM VciSci / ViSi, RDM BciSci / BiSi; Andrew 2026-11-12, was LBA
## VciSci / ViSi); columns LNR, RDM, RDM for the delta task's deltas, then the same for the CAF task's CAFs. Delta y axis on the left, CAF y axis on the right;
## y limits as plot_models6.R (delta -.0025 to .085, CAF 70-100); variant labels only, in Palatino (fonts embedded).
## Optional ROW3 adds a third row of one-mechanism models (Andrew, 2026-11-13): "LBA" = LNR Mi, LBA Vi, LBA Bi;
## "RDM" = LNR Mi, RDM Vi, RDM Bi. A row-3 panel whose architecture differs from its column's is labelled with it.
##   Rscript code/07_figures/plot_models6_mix.R [DELTA_TASK CAF_TASK [ROW3]]   (default Stroop ANT, no row 3)
## Needs derived/pp_raw/<label>_subj.rds (06_predict/pp_subj.R raw); writes figures/models6_<DELTA>_delta_<CAF>_caf.pdf
## (with ROW3: models6_<DELTA>_delta_<CAF>_caf_3row<ROW3>.pdf).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
a <- commandArgs(trailingOnly = TRUE)
TK <- if (length(a) >= 2) c(delta = a[1], caf = a[2]) else c(delta = "Stroop", caf = "ANT")
MODELS <- rbind(c("LNR_MciSci", "RDM_VciSci", "RDM_BciSci"),             # row 1: four mechanisms
                c("LNR_MiSi",   "RDM_ViSi",   "RDM_BiSi"))               # row 2: two mechanisms
ROW3 <- if (length(a) >= 3) toupper(a[3]) else ""
if (nzchar(ROW3)) { stopifnot(ROW3 %in% c("LBA", "RDM"))
  MODELS <- rbind(MODELS, c("LNR_Mi", paste0(ROW3, c("_Vi", "_Bi")))) }                # row 3: one mechanism
ARCHS <- sub("_.*", "", MODELS[1, ])
A <- read_tab("adequacy_display.csv")
YLIM <- list(delta = c(-.0025, .085), caf = c(70, 100))
FONT <- "Palatino"

OUT <- file.path(FIGDIR, emc2_name("models6", sprintf("%s_delta_%s_caf%s", TK[["delta"]], TK[["caf"]],
                                                     if (nzchar(ROW3)) paste0("_3row", ROW3) else "")))
pdf(OUT, width = 16, height = 6.6 * nrow(MODELS) / 2)
par(mfrow = c(nrow(MODELS), 6), oma = c(2.2, 2.8, 3.4, 2.8), mgp = c(1.4, 0.45, 0), tcl = -0.25)
for (meas in c("delta", "caf")) for (j in 1:3) {
  dat <- NULL
  for (i in 1:nrow(MODELS)) {
    r <- A[A$spec == paste(TK[[meas]], MODELS[i, j], sep = "_"), ][1, ]
    par(mfg = c(i, j + if (meas == "caf") 3 else 0),
        mar = if (meas == "delta") c(2.4, 2.4, 0.6, 0.4) else c(2.4, 0.4, 0.6, 2.4))
    if (!have_raw(r$label)) { blank("no predictives"); next }
    if (is.null(dat)) dat <- dat_of(r$label)
    pp <- readRDS(raw_file(r$label))
    if (meas == "caf") { invisible(draw_emc2(dat, pp, meas, YLIM[[meas]], yaxt = "n")); axis(4) }
    else invisible(draw_emc2(dat, pp, meas, YLIM[[meas]]))
    lab <- if (sub("_.*", "", MODELS[i, j]) == ARCHS[j]) sub("^[A-Z]+_", "", MODELS[i, j]) else sub("_", " ", MODELS[i, j])
    op <- par(family = FONT); label_emc2(lab, meas, cex = 0.95 * 1.5); par(op)
    cat(sprintf("%-6s %-5s %s %s\n", TK[[meas]], meas, r$label, r$class))
    rm(pp); gc(verbose = FALSE)
  }
}
cols <- (1:6 - 0.5) / 6
mtext(rep(ARCHS, 2), side = 3, outer = TRUE, line = 0.2, cex = 1.1, font = 2, at = cols)
mtext(sprintf("%s delta", TK[["delta"]]), side = 3, outer = TRUE, line = 1.8, cex = 1.15, font = 2, at = 0.25)
mtext(sprintf("%s CAF", TK[["caf"]]), side = 3, outer = TRUE, line = 1.8, cex = 1.15, font = 2, at = 0.75)
mtext("RT (I-C) (sec)", side = 2, outer = TRUE, line = 1.2, cex = 0.9)
mtext("CAF (%)", side = 4, outer = TRUE, line = 1.2, cex = 0.9)
mtext(c("Average RT (seconds)", "Bin Centre (%)"), side = 1, outer = TRUE, line = 0.9, cex = 0.85, at = c(0.25, 0.75))
invisible(dev.off())
if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded")   # Palatino is not a base PDF font
cat("wrote", OUT, "\n")
