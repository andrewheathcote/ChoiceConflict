## plot_mech_null.R transposed over all three tasks (Andrew, 2026-11-13): 3 x 6, rows the model's own predictives and
## <mechanism>=0 for each of its two i mechanisms (06_predict/pp_null.R); columns the deltas of Simon, Stroop and ANT,
## then their CAFs. One heading per row (left) and column (top), a y axis only on the leftmost delta and CAF columns
## (limits shared across each group: delta from a first drawing pass, CAF 70-100), one x
## label per group of three columns. EMC2 plot_delta / plot_caf; headings as plot_best5_grid.R.
##   Rscript code/07_figures/plot_mech_null_grid.R [MODEL]   (default RDM_ViSi)
## Needs derived/pp_raw/<label>_subj.rds and derived/pp_null/<label>_<par>0.rds; writes figures/mech_null_<MODEL>_grid.pdf.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
a <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(a) >= 1) a[1] else "RDM_ViSi"
MECH <- local({ v <- sub("^[A-Z]+_", "", MODEL); regmatches(v, gregexpr("[A-Z]i", v))[[1]] })
PAR  <- c(Vi = "v_iMd", Bi = "B_iMd", Mi = "m_iMd", Si = "s_iMd")[MECH]
D <- read_tab("display_jm.csv")
pp_null_path <- function(label, par) file.path(DERIVED, "pp_null", sprintf("%s_%s0.rds", label, par))
ROWS <- c(list(list(head = MODEL, file = function(lb) raw_file(lb))),
          lapply(seq_along(MECH), function(k) list(head = paste0(MECH[k], "=0"), file = function(lb) pp_null_path(lb, PAR[[k]]))))
LAB <- setNames(D$label[match(paste(TASKS, MODEL, sep = "_"), D$spec)], TASKS)

panels <- function(YL = NULL) {                         # draws every panel; returns the y range per measure
  rng <- list(delta = NULL, caf = NULL)
  for (ti in seq_along(TASKS)) {
    lb <- LAB[[TASKS[ti]]]; dat <- dat_of(lb)
    for (i in seq_along(ROWS)) {
      f <- ROWS[[i]]$file(lb); pp <- if (file.exists(f)) readRDS(f)
      for (meas in c("delta", "caf")) {
        j <- ti + if (meas == "caf") 3 else 0
        par(mfg = c(i, j), mar = c(2.2, if (ti == 1) 2.4 else 0.4, 0.5, 0.4))
        if (is.null(pp)) { blank("no predictives"); next }
        rng[[meas]] <- range(rng[[meas]], draw_emc2(dat, pp, meas, YL[[meas]], yaxt = if (ti == 1) "s" else "n"))
        if (!is.null(YL)) cat(sprintf("%-6s %-5s %-20s %s\n", TASKS[ti], meas, lb, basename(f)))
      }
      rm(pp); gc(verbose = FALSE)
    }
  }
  rng
}
layout_par <- function() par(mfrow = c(length(ROWS), 6), oma = c(2.4, 3.6, 1.8, 0.4), mgp = c(1.3, 0.45, 0), tcl = -0.25)
pdf(NULL); layout_par(); R <- panels(); invisible(dev.off())                # pass 1: shared y limits
YL <- list(delta = R$delta, caf = c(70, 100))

OUT <- fig_path(sprintf("mech_null_%s_grid.pdf", MODEL))
pdf(OUT, width = 16.54, height = 2.2 + 2.6 * length(ROWS))
layout_par(); invisible(panels(YL))
mtext(rep(TASKS, 2), side = 3, outer = TRUE, line = 0.3, cex = 1.15, font = 2, at = (1:6 - 0.5) / 6)
mtext(sapply(ROWS, `[[`, "head"), side = 2, outer = TRUE, line = 2.2, cex = 1.15, font = 2,
      at = (rev(seq_along(ROWS)) - 0.5) / length(ROWS))
mtext("RT (I-C) (sec)", side = 2, outer = TRUE, line = 0.6, cex = 0.95)
par(mfg = c(ceiling(length(ROWS) / 2), 4), mar = c(2.2, 2.4, 0.5, 0.4)); mtext("CAF (%)", side = 2, line = 1.4, cex = 0.95)   # beside the CAF axis
mtext(c("Average RT (seconds)", "Bin Centre (%)"), side = 1, outer = TRUE, line = 1.0, cex = 0.9, at = c(0.25, 0.75))
invisible(dev.off())
if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded")
cat("wrote", OUT, "\n")
