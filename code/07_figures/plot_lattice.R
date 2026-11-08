## Convergence grids, one page per task.
## Every mechanism belongs to a parameter and sits on a contrast: c = congruency,
## i = iMd. A parameter is therefore in one of four states (absent, c, i, ci),
## and an ensemble is a grid over those states. WDM/DDM mechanisms are
## congruency-only, so there each parameter is simply present or absent and the
## grid pairs them (a,z) against (v,s).
##
## DISPLAY VERSION (2026-10-20): copy of plot_lattice_jm.R. Outcome per cell from
## derived/tables/display_jm.csv (03_classify): converged; ridge-only
## shown as converged (same blue, hatched); failed (not displayed as converged:
## not converged, ridge-only Rhat >= 1.5, subject-mode caveat); fit running
## (campaign DDMs in play); not fit. mech from derived/tables/mech_key.csv by spec;
## never from file names.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
OUT   <- fig_path("convergence_grids.pdf")
TASKS <- c("Simon","Stroop","ANT")

## ---------------------------------------------------------------------------
## FIGURE CAPTION -- edit this text directly.
## ---------------------------------------------------------------------------
CAPTION <- paste(
"Each parameter is in one of four states: absent, c (constant across accumulators),",
"i (only the accumulator corresponding to the irrelevant attribute), or ci (both).",
"Rows and columns give those states.",
"\nPARAMETERS.  WDM and DDM: B = threshold, Z = bias, V = rate, S = sd(V).",
"LBA and RDM: B = threshold, V = rate, S = sd(V).  LNR: M = mean, S = sd.",
"\nDESIGN RULE: all possible single and double mechanism models, and all 3 and 4 mechanism models",
"that have at most one mechanism with only one of the c or i type.",
"\nAs the WDM and DDM have only one accumulator there is no c/i differentiation.",
"Models with more than four mechanisms are kept in <task>/<arch>/largek/.")
## ---------------------------------------------------------------------------
COL   <- c(conv = "#2B6CB0", ridge = "#2B6CB0", fail = "#C53030", run = "#F4F6F9", none = "#E2E8F0")

## joint-move family, best current fit per model (figs/display_jm.csv)
a  <- read_tab("display_jm.csv")
MK <- read_tab("mech_key.csv")              # one mechanism key (03_classify/mech_key.R)
a$variant <- sub("^[^_]+_", "", a$spec)
a$mech <- MK$mech[match(a$spec, MK$spec)]
stopifnot(!anyNA(a$mech), all(nzchar(a$mech)))
a$st <- ifelse(a$in_play, "run", ifelse(!a$display_converged, "fail",
               ifelse(a$class == "ridge-only, shown as converged", "ridge", "conv")))
a <- a[!(a$arch %in% c("WDM","DDM") & grepl("V", a$mech)), ]   # drop sv cases
RUNNING <- any(a$st == "run")
a <- a[!(a$arch %in% c("RDM","LBA") & nchar(a$mech) > 4), ]    # k<=4 everywhere

pstate <- function(m, up, lo) { l <- strsplit(m, "")[[1]]
  paste0(if (up %in% l) "c" else "", if (lo %in% l) "i" else "") }
bstate <- function(m, x, y)   { l <- strsplit(m, "")[[1]]
  paste0(if (x %in% l) x else "", if (y %in% l) y else "") }

cell <- function(d, rl, cl, rlab, clab, main, sub = "") {
  plot(NA, xlim = c(.5, 4.5), ylim = c(4.5, .5), axes = FALSE, xlab = "", ylab = "",
       main = "", cex.main = .92, font.main = 1)
  mtext(main, side = 3, line = 1.95, cex = 1.02, font = 1)
  if (nzchar(sub)) mtext(sub, side = 3, line = 1.45, cex = .62, col = "grey35")
  axis(3, 1:4, clab, tick = FALSE, line = -0.85, cex.axis = 1.35)
  axis(2, 1:4, rlab, tick = FALSE, line = -.55, las = 1, cex.axis = 1.35)
  for (i in 1:4) for (j in 1:4) {
    h  <- if (i == 1 && j == 1) integer(0) else which(d$r == rl[i] & d$c == cl[j])
    stopifnot(length(h) <= 1)
    st <- if (!length(h)) "none" else d$st[h]
    rect(j-.42, i-.42, j+.42, i+.42, col = COL[st], border = "white", lwd = 1.1)
    if (st == "run") rect(j-.42, i-.42, j+.42, i+.42, col = NA, border = "grey55", lwd = .9, lty = 2)
    if (st == "ridge") ridgehatch(j-.42, i-.42, j+.42, i+.42)   # slanted lines, as the heat grids (4 Nov; was a white dot)
    if (st == "fail") { segments(j-.2, i-.2, j+.2, i+.2, col = "white", lwd = 1.8)
                        segments(j-.2, i+.2, j+.2, i-.2, col = "white", lwd = 1.8) }
  }
}
ridgehatch <- function(x1, y1, x2, y2) rect(x1, y1, x2, y2, density = 22, angle = 45, col = "grey15", lwd = 0.5, border = NA)
legend_panel <- function() {
  plot(NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, xlab = "", ylab = "")
  text(.06, .86, "Fitting outcome", adj = 0, cex = 1.9, font = 1)
  ## "fit running" only while a campaign fit is in play (none since 3 Nov)
  ks <- if (RUNNING) 1:5 else c(1:3, 5); ly <- c(.68, .53, .38, .23, .08)[seq_along(ks)]; names(ly) <- ks
  for (k in ks) {
    rect(.08, ly[[as.character(k)]]-.055, .26, ly[[as.character(k)]]+.055, col = COL[k], border = "white", lwd = 1.2)
    text(.32, ly[[as.character(k)]], c("converged","ridge-only (DIC \u00b1250)","failed","fit running","not fit")[k],
         adj = 0, cex = 1.55)
  }
  ridgehatch(.08, ly[["2"]]-.055, .26, ly[["2"]]+.055)
  ly <- c(ly, "4" = .23)
  segments(.115, ly[["3"]]-.036, .225, ly[["3"]]+.036, col = "white", lwd = 2.2)
  segments(.115, ly[["3"]]+.036, .225, ly[["3"]]-.036, col = "white", lwd = 2.2)
  if (RUNNING) rect(.08, ly[["4"]]-.055, .26, ly[["4"]]+.055, col = NA, border = "grey55", lwd = .9, lty = 2)
}
SV <- c("","c","i","ci")
lab <- function(P) c("absent", paste0(P, c("c","i","ci")))

pdf(OUT, width = 13.5, height = 10)
for (t in TASKS) {
  par(mfrow = c(3,4), mar = c(1.2, 4.6, 4.2, 1), oma = c(6.5, 2, 5, 1))
  ## row 1 -- WDM, DDM, LNR, legend
  for (ar in c("WDM","DDM")) {
    s <- a[a$task == t & a$arch == ar, ]
    d <- data.frame(r = sapply(s$mech, bstate, "a","z"), c = sapply(s$mech, bstate, "v","s"), st = s$st)
    cell(d, c("","a","z","az"), c("","v","s","vs"),
         c("absent","B","Z","BZ"), c("absent","V","S","VS"), ar, "")
  }
  s <- a[a$task == t & a$arch == "LNR", ]
  d <- data.frame(r = sapply(s$mech, pstate, "M","m"), c = sapply(s$mech, pstate, "S","s"), st = s$st)
  cell(d, SV, SV, lab("M"), lab("S"), "LNR", "")
  legend_panel()
  ## rows 2 and 3 -- LBA then RDM, four noise-state panels each
  ## grid is V (rows) x S (columns); B varies across the four panels,
  ## running ci, c, i, absent -- the fullest specification first
  SP  <- c("ci", "c", "i", "")
  SPL <- c("ci", "c", "i", "absent")
  for (ar in c("LBA","RDM")) {
    s0 <- a[a$task == t & a$arch == ar, ]
    for (k in 1:4) {
      d <- data.frame(r = sapply(s0$mech, pstate, "V","v"), c = sapply(s0$mech, pstate, "S","s"),
                      p = sapply(s0$mech, pstate, "B","b"), st = s0$st)
      d <- d[d$p == SP[k], ]
      cell(d, SV, SV, lab("V"), lab("S"),
           if (k == 1) sprintf("%s:  B = %s", ar, SPL[k]) else sprintf("B = %s", SPL[k]), "")
    }
  }
  mtext(t, outer = TRUE, line = 2.6, cex = 1.5, font = 2)
  mtext(CAPTION, outer = TRUE, side = 1, line = 4.4, cex = .62, adj = 0, at = 0.012)
}
invisible(dev.off()); cat("wrote", OUT, "\n")
