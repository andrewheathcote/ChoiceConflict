## Choice recovery study scored by best-draw DIC, WAIC and PSIS-LOO (recch_loo.R outputs in loo/, DMC2 main f2f58985),
## with the real data's values (../../choice/compare/compare_choice.csv; its DIC column is standard DIC, so the best-draw
## DIC is rebuilt from meanD/Dmean/minD) and a DIC-grid-style figure: one panel per criterion, rows = data sets, cells =
## difference from the row's best (0 = best), true model boxed. No threshold is drawn (none is calibrated for choice).
##   Rscript code/08_recovery/choice/score_recch_ic.R    (writes recovery/choice/recch_ic_score.csv, figures/choice_recovery_ic.pdf)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
RC <- file.path(ROOT, "recovery", "choice")
FAM <- c("DDM", "LBA", "LNR", "RDM", "WDM"); CRIT <- c("DIC", "WAIC", "LOO")
bd_dic <- function(t) 2 * t$meanD - pmin(t$minD, t$Dmean)
L <- do.call(rbind, lapply(list.files(file.path(RC, "loo"), "[.]csv$", full.names = TRUE), read.csv, stringsAsFactors = FALSE))
dd <- read.csv(file.path(RC, "recch_dic.csv"), stringsAsFactors = FALSE); dd$set <- sub("_[A-Z]+(_v[01])?$", "", dd$spec)
L$DIC_dicjob <- dd$DIC[match(paste(L$set, L$model), paste(dd$set, dd$family))]
cat(sprintf("DIC in loo/ vs recch_dic.csv: max |diff| %.2f over %d fits\n", max(abs(L$DIC - L$DIC_dicjob)), nrow(L)))
Rl <- read_tab("compare_choice.csv")
Rl$DIC <- bd_dic(Rl)
## DIC for every set from recch_dic.csv (it exists before the slower WAIC/LOO jobs finish); WAIC/LOO from loo/
dd$model <- dd$family; dd <- merge(dd[, c("set", "model", "DIC")], L[, c("set", "model", "WAIC", "LOO")], all.x = TRUE)
V <- rbind(data.frame(set = "real", model = Rl$model, Rl[, CRIT]), dd[, c("set", "model", CRIT)])
SETS <- c("real", sprintf("rc_%s_r%d", rep(c("DDM", "LBA"), each = 5), rep(1:5, 2)))
lab <- c(real = "Real data", setNames(sprintf("%s-generated r%d", rep(c("DDM", "LBA"), each = 5), rep(1:5, 2)), SETS[-1]))
true_of <- function(s) if (s == "real") NA else sub("^rc_([A-Z]+)_r.*", "\\1", s)
## difference from the row best, per criterion
D <- lapply(setNames(CRIT, CRIT), function(cr) t(sapply(SETS, function(s) { v <- V[V$set == s, ]
  x <- setNames(v[[cr]], v$model)[FAM]; if (anyNA(x)) return(setNames(rep(NA_real_, 5), FAM)); x - min(x) })))
## score
S <- do.call(rbind, lapply(CRIT, function(cr) do.call(rbind, lapply(SETS, function(s) { d <- D[[cr]][s, ]
  if (all(is.na(d))) return(NULL); o <- order(d)
  data.frame(criterion = cr, set = s, true = true_of(s), winner = FAM[o[1]], margin = round(d[o[2]], 1),
             DDM_minus_LBA = round(d["LBA"] - d["DDM"], 1)) }))))
rownames(S) <- NULL; print(S, row.names = FALSE)
for (cr in CRIT) for (g in c("DDM", "LBA")) { s <- S[S$criterion == cr & S$true %in% g, ]
  cat(sprintf("%-4s %s-generated: recovered %d/%d, DDM - LBA %.0f to %.0f;  real %.0f\n", cr, g, sum(s$winner == g), nrow(s),
              min(s$DDM_minus_LBA), max(s$DDM_minus_LBA), S$DDM_minus_LBA[S$criterion == cr & S$set == "real"])) }
write_out(S, file.path(RC, "recch_ic_score.csv"))
## ---- figure, in the DIC grids' palette and log10(diff + 1) scale ----
RAMP <- colorRampPalette(rev(c("#313695","#4575B4","#74ADD1","#ABD9E9","#E0F3F8",
                               "#FFFFBF","#FEE090","#FDAE61","#F46D43","#D73027","#A50026")))(101)
HI <- max(unlist(D), na.rm = TRUE); tf <- function(v) log10(v + 1)
col_of <- function(v) if (is.na(v)) "grey88" else RAMP[round(100 * tf(v) / tf(HI)) + 1]
fmt <- function(v) if (is.na(v)) "pending" else if (v == 0) "0" else formatC(v, format = "f", digits = 0, big.mark = ",")
out <- fig_path("choice_recovery_ic.pdf")
pdf(out, width = 13, height = 6.2)
par(mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0)); plot.new(); plot.window(c(0, 1), c(0, 1), xaxs = "i", yaxs = "i")
nr <- length(SETS); top <- 0.86; rh <- 0.062; y <- top - (seq_len(nr) - 0.5) * rh - c(0, rep(0.02, 5), rep(0.04, 5))
lx <- 0.005; px0 <- 0.135; pw <- 0.255; gap <- 0.012; cw <- pw / 5
text(0.5, 0.965, "Choice model recovery: difference from the best model in each data set (0 = best; true model boxed)", font = 2, cex = 1.15)
text(lx, y, lab[SETS], adj = 0, cex = 0.9, font = ifelse(SETS == "real", 2, 1))
for (k in seq_along(CRIT)) { x0 <- px0 + (k - 1) * (pw + gap)
  text(x0 + pw / 2, 0.915, CRIT[k], font = 2, cex = 1.05)
  text(x0 + (1:5 - 0.5) * cw, top + 0.018, FAM, cex = 0.85)
  for (i in seq_len(nr)) { s <- SETS[i]; tr <- true_of(s)
    for (j in 1:5) { v <- D[[CRIT[k]]][s, j]; xc <- x0 + (j - 0.5) * cw
      rect(xc - cw * 0.47, y[i] - rh * 0.44, xc + cw * 0.47, y[i] + rh * 0.44, col = col_of(v), border = NA)
      if (!is.na(tr) && FAM[j] == tr) rect(xc - cw * 0.47, y[i] - rh * 0.44, xc + cw * 0.47, y[i] + rh * 0.44, border = "black", lwd = 2.2)
      tc <- if (is.na(v)) "grey40" else if (tf(v) / tf(HI) > 0.8) "white" else "black"
      text(xc, y[i], fmt(v), cex = 0.78, col = tc, font = if (!is.na(v) && v == 0) 2 else 1) } } }
## colour bar
bx0 <- 0.30; bx1 <- 0.70; by <- 0.035; tk <- c(0, 10, 100, 1000, 5000)
for (i in 1:101) rect(bx0 + (i - 1) / 101 * (bx1 - bx0), by, bx0 + i / 101 * (bx1 - bx0), by + 0.022, col = RAMP[i], border = NA)
tk <- tk[tk <= HI]; text(bx0 + tf(tk) / tf(HI) * (bx1 - bx0), by - 0.014, formatC(tk, format = "d", big.mark = ","), cex = 0.75)
text(bx0 - 0.01, by + 0.011, "difference from best (log scale)", adj = 1, cex = 0.8)
text(0.5, 0.105, "DIC: best-draw plug-in (as all project DICs); WAIC, LOO: DMC2 main f2f58985; deviance scale (-2 x elpd)", cex = 0.7, col = "grey30")
invisible(dev.off()); cat("wrote", out, "\n")
