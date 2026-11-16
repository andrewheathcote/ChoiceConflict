## The paper's Figure 14 (a diagram, no data): the supported race account, Flanker example: four pairs of accumulators (left responds <, right
## responds >): "Simple Choice" pairs for < and > stimuli and "Conflicted Choice" pairs for incongruent >><>> and
## congruent <<<<< (target < in both, so the left accumulator is correct). Lightning bolts from "AC (EEA & St)" go to
## the target-mismatch accumulators in simple choice and, in conflicted choice, the incongruent conflict-match and the
## congruent conflict-mismatch accumulators; bolts from "AC (Vi & Si)" (conflicted choice only) go to the
## conflict-mismatch accumulators.
##   Rscript code/07_figures/plot_attention.R      -> figures/attention_control.pdf
## Each accumulator shows one noisy accumulation path (schematic values): target match has the higher rate and target
## mismatch the larger noise. The conflicted-choice paths are the < simple-choice paths with the slope raised by D for
## the conflict-match and lowered by D for the conflict-mismatch accumulator (the direction of Vi > 0), so incongruent
## is slower and closer than simple choice and congruent faster. Paths stop when the first reaches threshold.
## Accumulator style follows the paper's Figures 3 and 6 (axis, threshold line, dotted baseline).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
out <- fig_path("attention_control.pdf")

W <- 1; H <- 1.6; GAP_ACC <- 0.6; GAP_PAIR <- 1.1; GAP_HALF <- 2.2           # accumulator size and spacing
pair_w <- 2 * W + GAP_ACC; half_w <- 2 * pair_w + GAP_PAIR
x0 <- c(0, pair_w + GAP_PAIR, half_w + GAP_HALF, half_w + GAP_HALF + pair_w + GAP_PAIR)   # left edge of each pair
acc_x <- function(p, k) x0[p] + (k - 1) * (W + GAP_ACC)                    # left edge of accumulator k (1 = <, 2 = >)
BOLT <- c("#F2B705", "#F07C24")                                           # control 1, 2

## schematic rates and noise for target match / mismatch, and the conflicted-choice slope change (conflict match +D,
## conflict mismatch -D)
V <- c(match = 2.3, mismatch = 1.2); S <- c(match = 0.3, mismatch = 0.6); D <- 0.3
## pair p: labels, stimulus, trial type, which accumulator matches the target (tm) and the conflict attribute (cm)
P <- list(list(lab = c("target\nmatch", "target\nmismatch"), stim = "<", type = "", tm = 1, cm = NA),
          list(lab = c("target\nmismatch", "target\nmatch"), stim = ">", type = "", tm = 2, cm = NA),
          list(lab = c("conflict\nmismatch", "conflict\nmatch"), stim = ">><>>", type = "", tm = 1, cm = 2),
          list(lab = c("conflict\nmatch", "conflict\nmismatch"), stim = "<<<<<", type = "", tm = 1, cm = 1))
## three levels above the accumulators. Simple choice: the figure header, "Simple Choice", AC (EEA & St). Conflicted
## choice: "Conflicted Choice", AC (Vi & Si), AC (EEA & St), so AC (EEA & St) is on the same level in both halves. In
## conflicted choice AC (Vi & Si) goes to the outer (conflict-mismatch) accumulators and AC (EEA & St) to the inner
## (conflict-match) one in the incongruent pair and the conflict-mismatch one in the congruent pair, which both reach;
## there the higher control's tip lands on the outer side (a third element of `to` shifts a tip), so no bolts cross. Each control has its centre x, level y, bolt colour, how far
## from its centre a bolt may start, and the accumulators acted on as c(pair, accumulator).
half_c <- function(h) (x0[2 * h - 1] + x0[2 * h] + pair_w) / 2
HEAD_Y <- 5.35; LEV2_Y <- 4.6; LEV3_Y <- 3.75
CTRL <- list(list(lab = "AC (EEA & St)", x = half_c(1), y = LEV3_Y, col = 1, reach = 0.8, to = list(c(1, 2), c(2, 1))),
             list(lab = "AC (Vi & Si)", x = half_c(2), y = LEV2_Y, col = 2, reach = 0.9, to = list(c(3, 1), c(4, 2, 0.2))),
             list(lab = "AC (EEA & St)", x = half_c(2), y = LEV3_Y, col = 1, reach = 0.6, to = list(c(3, 2), c(4, 2, -0.2))))

## one race: full-width paths for both accumulators on a dt grid (stop() then cuts them where the first reaches H)
DT <- 0.004; NT <- ceiling(W / DT)
paths <- function(p) { v <- s <- numeric(2)
  for (k in 1:2) { m <- if (k == P[[p]]$tm) "match" else "mismatch"; v[k] <- V[[m]]; s[k] <- S[[m]] }
  sapply(1:2, function(k) cumsum(c(0, rnorm(NT, v[k] * DT, s[k] * sqrt(DT))))) }
stop_at <- function(X) {
  hit <- apply(X, 2, function(x) which(x >= H)[1]); if (all(is.na(hit))) return(list(X = X, win = 0, t = Inf))
  end <- min(hit, na.rm = TRUE)
  list(X = X[seq_len(end), , drop = FALSE], win = which.min(replace(hit, is.na(hit), Inf)), t = (end - 1) * DT) }
## conflicted pair p from the < simple-choice paths X: conflict match slope +D, conflict mismatch -D
conflicted <- function(X, p) { tt <- (seq_len(nrow(X)) - 1) * DT
  for (k in 1:2) X[, k] <- X[, k] + (if (k == P[[p]]$cm) D else -D) * tt; stop_at(X) }
## the target-match accumulator wins, finishing within [t_lo, 0.9] of the width, with the loser visibly behind but in
## the race (between 0.3 and 0.85 of threshold) and never far below the baseline
ok <- function(r, p, t_lo = 0.55) r$win == P[[p]]$tm && r$t > t_lo && r$t < 0.9 && min(r$X) > -0.12 &&
  (function(x) x > 0.3 * H && x < 0.85 * H)(r$X[nrow(r$X), 3 - r$win])
## pair 2: first seed that is ok; pair 1: first seed ok for it and for both conflicted pairs made from it
pick <- function(p, also = function(X) TRUE) { for (seed in 1:20000) { set.seed(seed); X <- paths(p)
  if (ok(stop_at(X), p) && also(X)) { cat("pair", p, "seed", seed, "\n"); return(X) } }; stop("no seed for pair ", p) }
X1 <- pick(1, function(X) ok(conflicted(X, 3), 3, 0.4) && ok(conflicted(X, 4), 4, 0.4)); X2 <- pick(2)
R <- list(stop_at(X1), stop_at(X2), conflicted(X1, 3), conflicted(X1, 4))

accumulator <- function(x, path) {
  segments(x, 0, x, H); segments(x, H, x + W, H)                           # axis and threshold
  segments(x, 0, x + W, 0, lty = 3)                                         # baseline
  lines(x + (seq_along(path) - 1) * DT, path, lwd = 0.8)
}
## filled lightning bolt from (x1, y1) to its tip at (x2, y2); shape in (across, along) units, along 0..1
bolt <- function(x1, y1, x2, y2, col, w = 0.75) {
  s <- rbind(c(-0.06, 0), c(0.10, 0), c(0.02, 0.38), c(0.12, 0.38), c(-0.04, 1), c(0, 0.55), c(-0.10, 0.55))
  d <- c(x2 - x1, y2 - y1); L <- sqrt(sum(d^2)); u <- d / L; n <- c(-u[2], u[1])
  xy <- t(apply(s, 1, function(p) c(x1, y1) + p[2] * L * u + p[1] * w * n))
  polygon(xy, col = col, border = "grey25", lwd = 0.8)
}

pdf(out, width = 7.5, height = 3.65, family = "Helvetica")
par(mar = c(0.2, 0.2, 0.2, 0.2))
xr <- c(-0.2, max(x0) + pair_w + 0.2); plot(NA, xlim = xr, ylim = c(-0.7, 5.6), asp = 1, axes = FALSE, ann = FALSE)
LAB_Y <- H + 0.42
for (p in seq_along(P)) {
  for (k in 1:2) { accumulator(acc_x(p, k), R[[p]]$X[, k])
                   text(acc_x(p, k) + W / 2, LAB_Y, P[[p]]$lab[k], cex = 0.72) }
  xm <- x0[p] + pair_w / 2
  text(xm, -0.42, P[[p]]$stim, family = "mono", font = 2, cex = 1.25)
  if (nzchar(P[[p]]$type)) text(xm, -0.92, P[[p]]$type, cex = 0.72, font = 3)
}
text(x0[1], HEAD_Y, "Varieties of Attention Control (AC)", cex = 1.4, font = 2, adj = 0)
text(half_c(1), LEV2_Y, "Simple Choice", cex = 1.05, font = 2); text(half_c(2), HEAD_Y, "Conflicted Choice", cex = 1.05, font = 2)
for (c in CTRL) {
  text(c$x, c$y, c$lab, cex = 0.85, font = 2)
  for (t in c$to) {
    xt <- acc_x(t[1], t[2]) + W / 2 + if (length(t) > 2) t[3] else 0
    bolt(c$x + max(-c$reach, min(c$reach, 0.45 * (xt - c$x))), c$y - 0.22, xt, LAB_Y + 0.4, BOLT[c$col])
  }
}
dev.off(); cat("wrote", out, "\n")
