## What each mechanism does (Andrew, 2026-11-13): a model's own subject-level predictives beside predictives with one
## mechanism removed (06_predict/pp_null.R: its effect fixed at 0 for every subject and draw, same draws and seed).
## 2 x 3: deltas of the delta task (row 1) and CAFs of the CAF task (row 2); columns the fit (header the model),
## then <mechanism>=0 for each of its two i mechanisms. Drawn as plot_models6_mix.R (EMC2 plot_delta / plot_caf,
## Palatino); delta y limits shared across the row, from a first drawing pass; CAF 70-100 as plot_models6_mix.R. Column headers only (no panel labels).
##   Rscript code/07_figures/plot_mech_null.R [MODEL [DELTA_TASK CAF_TASK]]   (default RDM_ViSi Stroop ANT)
## Needs derived/pp_raw/<label>_subj.rds and derived/pp_null/<label>_<par>0.rds; writes figures/mech_null_<MODEL>.pdf
## (mech_null_<MODEL>_<DELTA>_delta_<CAF>_caf.pdf for other tasks).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
source(file.path(CODE, "lib", "emc2_plot_lib.R"))
a <- commandArgs(trailingOnly = TRUE)
MODEL <- if (length(a) >= 1) a[1] else "RDM_ViSi"
TK <- if (length(a) >= 3) c(delta = a[2], caf = a[3]) else c(delta = "Stroop", caf = "ANT")
## the model's two i mechanisms and their sampled parameters (mechanism letter -> parameter, as plot_effects.R)
MECH <- local({ v <- sub("^[A-Z]+_", "", MODEL); regmatches(v, gregexpr("[A-Z]i", v))[[1]] })
PAR  <- c(Vi = "v_iMd", Bi = "B_iMd", Mi = "m_iMd", Si = "s_iMd")[MECH]
D <- read_tab("display_jm.csv")
FONT <- "Palatino"
pp_null_path <- function(label, par) file.path(DERIVED, "pp_null", sprintf("%s_%s0.rds", label, par))
cols <- c(list(list(head = MODEL, file = function(lb) raw_file(lb))),
          lapply(seq_along(MECH), function(k) list(head = paste0(MECH[k], "=0"), file = function(lb) pp_null_path(lb, PAR[[k]]))))

panels <- function(YL = NULL) {                         # draws every panel; returns the y range per row
  rng <- list(delta = NULL, caf = NULL)
  for (meas in c("delta", "caf")) {
    lb <- D$label[D$spec == paste(TK[[meas]], MODEL, sep = "_")]; dat <- dat_of(lb); i <- if (meas == "delta") 1 else 2
    for (j in seq_along(cols)) {
      par(mfg = c(i, j))
      f <- cols[[j]]$file(lb); if (!file.exists(f)) { blank("no predictives"); next }
      pp <- readRDS(f)
      u <- draw_emc2(dat, pp, meas, YL[[meas]], yaxt = if (j > 1) "n" else "s")
      rng[[meas]] <- range(rng[[meas]], u)
      if (!is.null(YL)) { mtext(if (meas == "delta") "Average RT (seconds)" else "Bin Centre (%)", side = 1, line = 1.5, cex = 0.7)
        cat(sprintf("%-6s %-5s %-24s %s\n", TK[[meas]], meas, lb, basename(f))) }
      rm(pp); gc(verbose = FALSE)
    }
  }
  rng
}
layout_par <- function() par(mfrow = c(2, 3), oma = c(1, 4.2, 2.6, 0.6), mar = c(2.6, 2.4, 0.6, 0.4), mgp = c(1.4, 0.45, 0), tcl = -0.25)
pdf(NULL); layout_par(); YL <- panels(); invisible(dev.off())              # pass 1: shared delta y limits
YL$caf <- c(70, 100)                                                       # CAF as plot_models6_mix.R

OUT <- fig_path(sprintf("mech_null_%s%s.pdf", MODEL,
                        if (identical(unname(TK), c("Stroop", "ANT"))) "" else sprintf("_%s_delta_%s_caf", TK[["delta"]], TK[["caf"]])))
pdf(OUT, width = 9.5, height = 6.4)
layout_par(); invisible(panels(YL))
mtext(sapply(cols, `[[`, "head"), side = 3, outer = TRUE, line = 0.6, cex = 1.15, font = 2, at = (seq_along(cols) - 0.5) / 3, family = FONT)
mtext(c(sprintf("%s CAF (%%)", TK[["caf"]]), sprintf("%s delta: RT (I-C) (sec)", TK[["delta"]])), side = 2, outer = TRUE,
      line = 1.4, cex = 0.9, at = c(0.27, 0.77))
invisible(dev.off())
if (nzchar(Sys.which("gs"))) embedFonts(OUT) else warning("Ghostscript not found: fonts not embedded")
cat("wrote", OUT, "\n")
