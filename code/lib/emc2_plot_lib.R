## Shared by the *_emc2.R figure scripts: the display figures redrawn with EMC2's own plot_delta / plot_caf, as the
## archive's plot_best5.R / plot_figA.R / plot_figB.R drew the legacy figures (same arguments: delta on correct RTs
## at the 10-90% percentiles, CAF smoothed over a 5% window). EMC2 pools each draw's trials over subjects, so these
## curves are not the subject means of the *_display.R figures.
## Predictives: derived/pp_raw/<label>_subj.rds (06_predict/pp_subj.R raw), subject level only.
## Drawn with the predictive build (config.R BUILD["predict"]), so the figures do not depend on a default library.
use_build(BUILD[["predict"]])
source(file.path(CODE, "lib", "pp_display_lib.R"))
PCT   <- c(10, 30, 50, 70, 90)                     # the five quantiles the misfit statistics use
CIfun <- function(d) factor(d$RS != d$IS, labels = c("Congruent", "Incongruent"))
## accuracy-coded data (Stroop WDM/DDM, dat_acc.RData) carries R as FALSE/TRUE
mkacc  <- function(acc) if (acc) function(d) d$R == "TRUE" else function(d) d$RS == d$R
mkcfun <- function(acc) { ok <- mkacc(acc)
  function(d) factor(!ok(d), levels = c(FALSE, TRUE), labels = c("Correct", "Error")) }
FM <- read_tab("FITS.csv")
dat_of <- function(label) spec_dat(read_spec(FM$spec[FM$label == label]))   # the data the fit was made on, via its spec
raw_file <- function(label) pp_raw_path(label)
have_raw <- function(label) !is.null(label) && !is.na(label) && file.exists(raw_file(label))
emc2_suf <- paste0("_", LEVEL)
## figure file name, the simplest (Andrew, 2026-10-21): best5_Simon.pdf at subject level (the paper's delta / CAF
## figures), best5_Simon_hyper.pdf at hyper level
emc2_name <- function(stem, rest) sprintf("%s_%s%s.pdf", stem, rest, if (LEVEL == "hyper") "_hyper" else "")

## one panel; returns par("usr")[3:4] so callers can share y limits (the archive scripts' two-pass scheme);
## ... goes to plot_delta / plot_caf and on to plot() (e.g. yaxt = "n" to draw the y axis elsewhere)
draw_emc2 <- function(dat, pp, meas, YL = NULL, ...) {
  acc <- identical(levels(factor(dat$R)), c("FALSE", "TRUE")); yl <- c(if (is.null(YL)) list() else list(ylim = YL), list(...))
  if (meas == "delta")
    do.call(plot_delta, c(list(dat, pp, delta_factor = "CI", main = "", layout = NULL, legendpos = NA,
      functions = list(CI = CIfun, C = mkcfun(acc)), add_percentiles = PCT, xlab = "", ylab = ""), yl))
  else
    do.call(plot_caf, c(list(dat, pp, caf_factor = "CI", main = "", layout = NULL, functions = list(CI = CIfun),
      smooth_window = 5, accuracy_function = mkacc(acc), legendpos = c(NA, NA), xlab = "", ylab = ""), yl))
  par("usr")[3:4]
}
## panel label: top left on delta, bottom right on CAF (the archive scripts' placement); bottom = TRUE puts a delta
## label bottom left, clear of a curve that reaches the top of its panel
label_emc2 <- function(txt, meas, bottom = FALSE, cex = 0.95) { u <- par("usr")
  if (meas == "delta" && bottom) text(u[1] + 0.04 * (u[2] - u[1]), u[3] + 0.07 * (u[4] - u[3]), txt, adj = 0, cex = cex, font = 2)
  else if (meas == "delta") text(u[1] + 0.04 * (u[2] - u[1]), u[4] - 0.06 * (u[4] - u[3]), txt, adj = 0, cex = cex, font = 2)
  else text(u[2] - 0.04 * (u[2] - u[1]), u[3] + 0.07 * (u[4] - u[3]), txt, adj = 1, cex = cex, font = 2) }
