## Shared code for the goodness-of-fit figures on the display fits (plot_best5_display.R,
## plot_models6_display.R, plot_split6_display.R). Reads the posterior-predictive primitives in
## derived/pp_subj/<label>.rds (06_predict/pp_subj.R: 100 draws x 518 subjects).
## The group curve is an aggregate over subjects of each subject's five correct-RT quantiles
## (10-90%) by congruency and five-bin accuracy. AGG chooses the aggregate for the RT quantiles:
##   mean    (default) the Vincentised mean, identical to the stored gqc/gqi
##   median  robust to a few runaway simulated subjects (HANDOFF.md "Known weakness")
## Set it with the environment variable PP_AGG=mean|median. Accuracy is always the mean over
## subjects: it is bounded, so it has no runaway problem, and a subject's bin accuracy rests on a
## few trials, so it is discrete and mostly 1 -- the median over subjects is then 1 for data and
## model alike and the CAF misfit collapses to 0 (found 2026-10-20).
## Subject-level predictives only (Andrew, 4 Nov: a project rule; hyper versions removed): pp_subj/<label>.rds
## (pp_subj.R), each real subject simulated from its own alpha draws. PP_LEVEL is kept only so that old calls
## with PP_LEVEL=subj still work; any other value stops.
LEVEL <- Sys.getenv("PP_LEVEL", unset = "subj")
if (LEVEL != "subj") stop("subject-level posterior predictives only in this project (Andrew, 4 Nov); hyper is not made")
AGG <- Sys.getenv("PP_AGG", unset = "mean"); stopifnot(AGG %in% c("mean", "median"))
meanf <- function(m) colMeans(m, na.rm = TRUE)
aggf  <- if (AGG == "mean") meanf else function(m) apply(m, 2, median, na.rm = TRUE)
PCT  <- c(10, 30, 50, 70, 90)
COLC <- "#2a78d6"; COLI <- "#eb6834"; COLD <- "grey35"     # congruent, incongruent, delta (validated pair)
band <- function(col, a = 0.22) adjustcolor(col, alpha.f = a)
HALF <- function(v) mean(v[4:5]) - mean(v[1:2])            # slope rule of the archive's rank_subjects.R

read_summ <- function(dest) readRDS(file.path(DERIVED, "pp_subj", dest))
have_pp <- function(dest) !is.na(dest) && file.exists(file.path(DERIVED, "pp_subj", dest))
## group curves from per-subject arrays; `sel` (draw x subject logical, or subject logical for the
## data) restricts the subjects, as the split figure needs
grp_pp <- function(pp, sel = NULL) {
  nd <- dim(pp$sqc)[1]; ok <- !is.na(pp$gqc[, 1])
  one <- function(a, k, fn) { m <- a[k, , ]; if (!is.null(sel)) m <- m[if (is.matrix(sel)) sel[k, ] else sel, , drop = FALSE]
    fn(m) }
  f <- function(a, fn) t(vapply(which(ok), function(k) one(a, k, fn), numeric(dim(a)[3])))
  list(qc = f(pp$sqc, aggf), qi = f(pp$sqi, aggf), aI = f(pp$scafI, meanf), aC = f(pp$scafC, meanf))
}
grp_obs <- function(ob, sel = NULL) {
  g <- function(m, fn) { if (!is.null(sel)) m <- m[sel, , drop = FALSE]; fn(m) }
  list(qc = g(ob$sqc, aggf), qi = g(ob$sqi, aggf), aI = g(ob$scafI, meanf), aC = g(ob$scafC, meanf))
}
q3 <- function(m) apply(m, 2, quantile, c(.025, .5, .975), na.rm = TRUE)
curves_from <- function(P, O) list(
  dx = apply((P$qc + P$qi) / 2, 2, median), dy = q3(P$qi - P$qc),
  ox = (O$qc + O$qi) / 2, oy = O$qi - O$qc,
  cC = 100 * q3(P$aC), cI = 100 * q3(P$aI), oC = 100 * O$aC, oI = 100 * O$aI)
curves <- function(dest) { z <- read_summ(dest); curves_from(grp_pp(z$pp), grp_obs(z$obs)) }

blank <- function(msg) { plot(NA, xlim = 0:1, ylim = 0:1, axes = FALSE, xlab = "", ylab = "")
  box(col = "grey85"); text(.5, .55, msg, col = "grey45", cex = .95) }
## class tag for a panel label
tag <- function(cl) switch(cl, "converged" = "", "ridge-only, shown as converged" = " (ridge-only)",
                           " (not converged)")
draw_delta <- function(cv, xl, yl, label, legend = FALSE) {
  plot(NA, xlim = xl, ylim = yl, xlab = "", ylab = "", las = 1)
  abline(h = 0, col = "grey80", lwd = .8)
  polygon(c(cv$dx, rev(cv$dx)), c(cv$dy[1, ], rev(cv$dy[3, ])), col = band(COLD), border = NA)
  lines(cv$dx, cv$dy[2, ], col = COLD, lwd = 2)
  lines(cv$ox, cv$oy, col = "black", lwd = .8, lty = 3)
  points(cv$ox, cv$oy, pch = 21, bg = "black", col = "white", cex = 1.25, lwd = .8)
  u <- par("usr"); text(u[1] + .04 * diff(u[1:2]), u[4] - .06 * diff(u[3:4]), label, adj = 0, cex = .95, font = 2)
  if (legend) legend("bottomright", bty = "n", cex = .85, inset = c(.01, .02),
    legend = c("data", "model median", "model 95%"), pch = c(21, NA, 15), pt.bg = "black",
    col = c("white", COLD, band(COLD, .4)), lwd = c(NA, 2, NA), pt.cex = c(1.2, 1, 2.2))
}
draw_caf <- function(cv, yl, label, legend = FALSE) {
  plot(NA, xlim = c(0, 100), ylim = yl, xlab = "", ylab = "", las = 1, xaxt = "n")
  axis(1, at = PCT)
  for (k in list(list(p = cv$cC, o = cv$oC, col = COLC), list(p = cv$cI, o = cv$oI, col = COLI))) {
    polygon(c(PCT, rev(PCT)), c(k$p[1, ], rev(k$p[3, ])), col = band(k$col), border = NA)
    lines(PCT, k$p[2, ], col = k$col, lwd = 2)
    lines(PCT, k$o, col = k$col, lwd = .8, lty = 3)
    points(PCT, k$o, pch = 21, bg = k$col, col = "white", cex = 1.25, lwd = .8)
  }
  u <- par("usr"); text(u[2] - .04 * diff(u[1:2]), u[3] + .07 * diff(u[3:4]), label, adj = 1, cex = .95, font = 2)
  if (legend) legend("bottomleft", bty = "n", cex = .85, inset = c(.01, .14),
    legend = c("congruent", "incongruent"), pch = 16, col = c(COLC, COLI),
    lwd = 2, lty = 1, pt.cex = 1.1, seg.len = 1.6)
}
## shared delta limits over a list of curve sets
dlims <- function(cvs) {
  cvs <- Filter(Negate(is.null), cvs)
  yl <- range(unlist(lapply(cvs, function(c) c(c$dy, c$oy))), na.rm = TRUE)
  xl <- range(unlist(lapply(cvs, function(c) c(c$dx, c$ox))), na.rm = TRUE)
  list(x = xl + c(-.04, .04) * diff(xl), y = yl + c(-.06, .10) * diff(yl))
}
## output file suffix: nothing for the mean (the default), "_median" otherwise
aggsuf <- if (AGG == "mean") "" else "_median"
levsuf <- paste0("_", LEVEL, aggsuf)                         # figure file suffix, e.g. _hyper, _subj_median
## figure file name: one set of PDFs with the simplest names (Andrew, 2026-10-21). The paper's delta / CAF figures
## are the EMC2 versions (*_emc2.R); running these Vincentised scripts overwrites them.
figname <- function(stem, rest) sprintf("%s_%s.pdf", stem, rest)
