## Score the choice recovery study by best-draw DIC and BPIC (recch_dic.csv from dic.R): per simulated data set and
## criterion, the winner among the five candidates, its margin over the runner-up, the true model's gap to the winner and
## the DDM - LBA difference (LBA minus DDM; positive = DDM ahead), set against the real data's values for the same
## candidates (derived/tables/choice_dic.csv, the fits in CHOICE_FULL). Also a DIC-grid-style figure: one panel per
## criterion, rows = data sets, cells = difference from the row's best (0 = best), true model boxed.
##   Rscript code/08_recovery/choice/score.R     -> recovery/choice/recch_score.csv, figures/choice_recovery.pdf
## From newsampler/gadi/choice_rec/score_recch.R and score_recch_ic.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
RC <- file.path(ROOT, "recovery", "choice")
FAM <- c("DDM", "LBA", "LNR", "RDM", "WDM"); CRIT <- c("DIC", "BPIC")
D <- read.csv(file.path(RC, "recch_dic.csv"), stringsAsFactors = FALSE); D$set <- sub("_[A-Z]+(_v[01])?$", "", D$spec)
real <- read_tab("choice_dic.csv"); real <- real[match(CHOICE_FULL[FAM], real$label), ]; stopifnot(!anyNA(real$label))
V <- rbind(data.frame(set = "real", model = FAM, real[, CRIT]), data.frame(set = D$set, model = D$family, D[, CRIT]))
SETS <- c("real", sprintf("rc_%s_r%d", rep(c("DDM", "LBA"), each = 5), rep(1:5, 2)))
lab <- c(real = "Real data", setNames(sprintf("%s-generated r%d", rep(c("DDM", "LBA"), each = 5), rep(1:5, 2)), SETS[-1]))
true_of <- function(s) if (s == "real") NA else sub("^rc_([A-Z]+)_r.*", "\\1", s)
unmet <- function(s) paste(D$family[D$set == s & D$stop_rule_met != "TRUE"], collapse = ",")
## difference from the row best, per criterion
G <- lapply(setNames(CRIT, CRIT), function(cr) t(sapply(SETS, function(s) { v <- V[V$set == s, ]
  x <- setNames(v[[cr]], v$model)[FAM]; stopifnot(!anyNA(x)); x - min(x) })))
S <- do.call(rbind, lapply(CRIT, function(cr) do.call(rbind, lapply(SETS, function(s) { d <- G[[cr]][s, ]; o <- order(d)
  data.frame(criterion = cr, set = s, generator = true_of(s), winner = FAM[o[1]], runner_up = FAM[o[2]],
             margin = round(d[o[2]], 1), true_gap = if (s == "real") NA else round(d[true_of(s)], 1),
             DDM_minus_LBA = round(d["LBA"] - d["DDM"], 1), order = paste(FAM[o], collapse = " < "),
             unmet = if (s == "real") "" else unmet(s)) }))))
rownames(S) <- NULL; print(S, row.names = FALSE)
for (cr in CRIT) for (g in c("DDM", "LBA")) { s <- S[S$criterion == cr & S$generator %in% g, ]
  cat(sprintf("%-4s %s-generated: recovered %d/%d, margin %.0f to %.0f, DDM - LBA %.0f to %.0f;  real %.0f\n", cr, g,
              sum(s$winner == g), nrow(s), min(s$margin), max(s$margin), min(s$DDM_minus_LBA), max(s$DDM_minus_LBA),
              S$DDM_minus_LBA[S$criterion == cr & S$set == "real"])) }
write_out(S, file.path(RC, "recch_score.csv"))
## ---- figure, in the DIC grids' palette and log10(diff + 1) scale ----
RAMP <- colorRampPalette(rev(c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
                               "#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")))(101)
HI <- max(unlist(G)); tf <- function(v) log10(v + 1)
col_of <- function(v) RAMP[round(100 * tf(v) / tf(HI)) + 1]
fmt <- function(v) if (v == 0) "0" else formatC(v, format = "f", digits = 0, big.mark = ",")
out <- fig_path("choice_recovery.pdf")
pdf(out, width = 9.5, height = 6.2)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
nr <- length(SETS); top <- 0.86; rh <- 0.062; y <- top - (seq_len(nr) - 0.5) * rh - c(0, rep(0.02, 5), rep(0.04, 5))
lx <- 0.008; px0 <- 0.185; pw <- 0.395; gap <- 0.015; cw <- pw / 5
text(0.5, 0.965, "Choice model recovery: difference from the best model in each data set (0 = best; true model boxed)", font = 2, cex = 1.05)
text(lx, y, lab[SETS], adj = 0, cex = 0.9, font = ifelse(SETS == "real", 2, 1))
for (k in seq_along(CRIT)) { x0 <- px0 + (k - 1) * (pw + gap)
  text(x0 + pw / 2, 0.915, CRIT[k], font = 2, cex = 1.05)
  text(x0 + (1:5 - 0.5) * cw, top + 0.018, FAM, cex = 0.85)
  for (i in seq_len(nr)) { s <- SETS[i]; tr <- true_of(s)
    for (j in 1:5) { v <- G[[CRIT[k]]][s, j]; xc <- x0 + (j - 0.5) * cw
      rect(xc - cw * 0.47, y[i] - rh * 0.44, xc + cw * 0.47, y[i] + rh * 0.44, col = col_of(v), border = NA)
      if (!is.na(tr) && FAM[j] == tr) rect(xc - cw * 0.47, y[i] - rh * 0.44, xc + cw * 0.47, y[i] + rh * 0.44, border = "black", lwd = 2.2)
      text(xc, y[i], fmt(v), cex = 0.78, col = if (tf(v) / tf(HI) > 0.8) "white" else "black", font = if (v == 0) 2 else 1) } } }
## colour bar
bx0 <- 0.36; bx1 <- 0.76; by <- 0.035; tk <- c(0, 10, 100, 1000, 5000)
for (i in 1:101) rect(bx0 + (i - 1) / 101 * (bx1 - bx0), by, bx0 + i / 101 * (bx1 - bx0), by + 0.022, col = RAMP[i], border = NA)
tk <- tk[tk <= HI]; text(bx0 + tf(tk) / tf(HI) * (bx1 - bx0), by - 0.014, formatC(tk, format = "d", big.mark = ","), cex = 0.75)
text(bx0 - 0.01, by + 0.011, "difference from best (log scale)", adj = 1, cex = 0.8)
text(0.5, 0.105, "Best-draw conditional DIC and BPIC (as all project DICs), EMC2 4f1fa123", cex = 0.7, col = "grey30")
invisible(dev.off()); cat("wrote", out, "\n")
