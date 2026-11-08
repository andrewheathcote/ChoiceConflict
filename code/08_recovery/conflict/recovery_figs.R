## Figures of the conflict recovery studies, from recovery/conflict/*.csv:
##   Rscript code/08_recovery/conflict/recovery_figs.R    -> figures/recovery_jm_fig.pdf, figures/recovery_rec5_fig.pdf
## recovery_jm_fig: left, the true-vs-winner DIC gap of every data set of studies 1-3 (LNR, RDM, LBA), with the DIC floor;
##   right, the cross-architecture study's margins: the true model's lead over the nearest same-architecture candidate
##   and over the best candidate of the other architecture.
## recovery_rec5_fig: per cell, DIC(best two-mechanism candidate) - DIC(one-mechanism model) on each simulated set from
##   the one-mechanism generator (G1) and the two-mechanism generator (G2), and on the real data. Negative favours two.
## From newsampler/recovery_jm_fig.R (its legacy-sampler points removed) and rec5/recovery_rec5_fig.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
RC <- file.path(ROOT, "recovery", "conflict")
FLOOR <- read.csv(file.path(INPUTS, "noise_floors.csv"))$floor[1]
s <- read.csv(file.path(RC, "recovery_jm_sets.csv"), stringsAsFactors = FALSE); m <- read.csv(file.path(RC, "recovery_jm.csv"), stringsAsFactors = FALSE)
pdf(fig_path("recovery_jm_fig.pdf"), width = 9, height = 4.2); par(mfrow = c(1, 2), mar = c(4, 4.2, 2, 1), las = 1, cex = 0.85)
ST <- c("LNR", "RDM", "LBA"); off <- c(LNR = 0, RDM = 22, LBA = 44)
plot(NA, xlim = c(0, 64), ylim = c(0, 1.2 * FLOOR), xaxt = "n", xlab = "", ylab = "DIC(true) - DIC(winner)", main = "Gap when the generator lost")
for (st in ST) { x <- s[s$study == st, ]; i <- seq_len(nrow(x)) + off[[st]]; points(i, x$gap, pch = 19, col = "#b2182b") }
axis(1, at = off + 10.5, labels = ST, tick = FALSE)
abline(h = FLOOR, lty = 2, col = "#b2182b"); text(64, FLOOR, paste("DIC floor", FLOOR), adj = c(1, -0.4), col = "#b2182b", cex = 0.85)
x <- m[m$study == "XA", ]; x$arch <- sub("_.*", "", x$candidate)
r <- do.call(rbind, lapply(split(x, x$data_file), function(z) { ta <- z$arch[z$true_model]; t <- z$DIC[z$true_model]
  data.frame(gen = z$generator[1], rep = z$replicate[1], within = min(z$DIC[!z$true_model & z$arch == ta]) - t,
             across = min(z$DIC[z$arch != ta]) - min(z$DIC[z$arch == ta])) }))
r <- r[order(r$gen, r$rep), ]; i <- seq_len(nrow(r))
plot(i, r$across, log = "y", ylim = c(20, 3000), pch = 19, col = "#2166ac", xaxt = "n", xlab = "", ylab = "true model's DIC lead",
     main = sprintf("Cross-architecture: %d/%d recovered", sum(s$recovered[s$study == "XA"], na.rm = TRUE), sum(s$study == "XA")))
points(i, r$within, pch = 17, col = "#b2182b"); abline(h = FLOOR, lty = 2, col = "#b2182b"); abline(v = 5.5, col = "grey70")
axis(1, at = c(3, 8), labels = c("DDM_BiVc data", "RDM_ViSi data"), tick = FALSE)
legend("left", c("over the other architecture's best", "over the nearest same-architecture model"), pch = c(19, 17),
       col = c("#2166ac", "#b2182b"), bty = "n", cex = 0.8)
invisible(dev.off())

s5 <- read.csv(file.path(RC, "recovery_rec5_sets.csv"), stringsAsFactors = FALSE); cl <- read.csv(file.path(RC, "recovery_rec5_cells.csv"), stringsAsFactors = FALSE)
cl$lab <- sprintf("%s %s: %s vs %s", cl$task, cl$arch, cl$G1, cl$G2); y <- rev(seq_len(nrow(cl))); names(y) <- cl$cell
pdf(fig_path("recovery_rec5_fig.pdf"), width = 8, height = 4.2); par(mar = c(4, 12, 1, 1), las = 1, cex = 0.85)
plot(NA, xlim = c(-450, 450), ylim = c(0.5, nrow(cl) + 0.5), yaxt = "n", ylab = "",
     xlab = "DIC(best two-mechanism) - DIC(one-mechanism)    (negative favours two mechanisms)")
abline(v = 0, col = "grey70"); abline(h = y, col = "grey92")
axis(2, at = y, labels = cl$lab, tick = FALSE)
g1 <- s5[s5$gen == "g1", ]; g2 <- s5[s5$gen == "g2", ]
points(g1$best2_minus_g1, y[g1$cell] + 0.15, pch = 1, col = "#2166ac")
points(g2$best2_minus_g1, y[g2$cell] - 0.15, pch = 2, col = "#b2182b")
points(cl$real_best2_minus_g1, y[cl$cell], pch = 18, cex = 1.7)
legend("bottomright", c("simulated from the one-mechanism model", "simulated from the two-mechanism model", "real data"),
       pch = c(1, 2, 18), col = c("#2166ac", "#b2182b", "black"), bty = "n", cex = 0.8)
invisible(dev.off()); cat("wrote", fig_path("recovery_jm_fig.pdf"), "and", fig_path("recovery_rec5_fig.pdf"), "\n")
