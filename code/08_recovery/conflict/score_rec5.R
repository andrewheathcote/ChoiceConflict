## Recovery study 5: one mechanism or two (PLAN_recovery_onemech.md). Scores the jmh-rec5 fits:
##   Rscript code/08_recovery/conflict/score_rec5.R   (reads <work>/summ/<spec>_jmhrec5.rds from 02_fit/summarise_fit.R run
##   with CC_DERIVED=<work>, the fit records <work>/fits/<spec>_jmhrec5.txt, and <work>/screen_{old,fixed}.rds when present)
## Per data set (cell x generator x replicate, 4 candidates): the conditional-DIC winner, whether it is the generator,
## the true-vs-winner gap, and d21 = DIC(G2) - DIC(G1) (negative: the two-mechanism model looks better). Per cell:
##   false alarms  - G1 data won by a two-mechanism candidate, and the largest such margin;
##   misses        - G2 data won by G1, and the largest such margin;
##   real gap      - the real data's d21 (figs/display_jm.csv) placed against the simulated d21 under each generator;
##   floor         - the largest true-vs-winner gap in the cell, for comparison with the DIC floor (150).
## RDM fits (c1, c5, c7) were made with the pre-26-Oct C++ RDM likelihood; a fit with a subject whose log-likelihood at
## the posterior mean is > 5 above the fixed build's is flagged (dic_gain = 2 x summed inflation, as compare_ll.R).
## Writes recovery/conflict/recovery_rec5.csv (per fit), recovery_rec5_sets.csv and recovery_rec5_cells.csv.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work"))
R5 <- file.path(ROOT, "recovery", "conflict"); FLOOR <- 150
map <- read.csv(file.path(WORK, "rec5_map.csv"), stringsAsFactors = FALSE)
CELLS <- readRDS(file.path(WORK, "cells.rds"))
rd <- function(f) { x <- readLines(f); k <- sub(":.*$", "", x); v <- sub("^[^:]*: ", "", x)
  keep <- !duplicated(k, fromLast = TRUE); setNames(v[keep], k[keep]) }
map$DIC <- NA_real_; map$rule <- NA_real_; map$met <- NA; map$meanD <- map$Dmean <- map$minD <- NA_real_
for (i in seq_len(nrow(map))) {
  fs <- file.path(WORK, "summ", paste0(map$spec[i], "_jmhrec5.rds")); ft <- file.path(WORK, "fits", paste0(map$spec[i], "_jmhrec5.txt"))
  if (file.exists(fs)) { s <- readRDS(fs); if (is.data.frame(s$dic)) { map$DIC[i] <- s$dic$DIC[1]
    map$meanD[i] <- s$dic$meanD[1]; map$Dmean[i] <- s$dic$Dmean[1]; map$minD[i] <- s$dic$minD[1] } }
  if (file.exists(ft)) { r <- rd(ft); g <- function(k) if (k %in% names(r)) r[[k]] else NA
    map$rule[i] <- suppressWarnings(as.numeric(g("stop_rule_value"))); map$met[i] <- g("stop_rule_met") %in% "TRUE" }
}
## likelihood screen (RDM only)
map$ll_n_up <- NA_integer_; map$ll_dic_gain <- NA_real_
fo <- file.path(WORK, "screen_old.rds"); fn <- file.path(WORK, "screen_fixed.rds")
if (file.exists(fo) && file.exists(fn)) { o <- readRDS(fo); n <- readRDS(fn)
  for (sn in intersect(names(o), names(n))) { if (is.null(o[[sn]]$ll) || is.null(n[[sn]]$ll)) next
    s <- intersect(names(o[[sn]]$ll), names(n[[sn]]$ll)); d <- o[[sn]]$ll[s] - n[[sn]]$ll[s]; up <- d > 5
    k <- map$spec == sn; map$ll_n_up[k] <- sum(up); map$ll_dic_gain[k] <- round(2 * sum(d[up])) } }
map$DIC_standard <- 2 * map$meanD - map$Dmean          # posterior-mean plug-in, for dic_variants.R
map$set <- sprintf("%s_%s_r%d", map$cell, map$gen, map$replicate)
write_out(map[order(map$spec), ], file.path(R5, "recovery_rec5.csv"))

## per data set
sets <- do.call(rbind, lapply(split(map, map$set), function(x) {
  cl <- CELLS[[x$cell[1]]]; dic <- setNames(x$DIC, x$candidate); done <- !anyNA(dic)
  w <- if (done) names(which.min(dic)) else NA
  data.frame(cell = x$cell[1], task = x$task[1], arch = x$arch[1], gen = x$gen[1], generator = x$generator[1], replicate = x$replicate[1],
             n_done = sum(!is.na(dic)), n_met = sum(x$met %in% TRUE), true_met = isTRUE(x$met[x$true_model]),
             winner = w, winner_mech = if (done) x$n_mech[x$candidate == w] else NA, recovered = if (done) x$true_model[x$candidate == w] else NA,
             gap = if (done) dic[x$true_model] - min(dic) else NA,
             d21 = dic[cl$g[2]] - dic[cl$g[1]],                                      # G2 - G1
             best2_minus_g1 = if (done) min(dic[x$n_mech == 2]) - dic[cl$g[1]] else NA,
             ll_flag = any(x$ll_n_up > 0, na.rm = TRUE), row.names = NULL)
}))
sets <- sets[order(sets$cell, sets$gen, sets$replicate), ]
write_out(sets, file.path(R5, "recovery_rec5_sets.csv"))

## per cell, with the real data's gaps
D <- read_tab("display_jm.csv")
realDIC <- function(task, arch, cand) { v <- D$DIC[D$spec == paste(task, arch, cand, sep = "_")]; if (length(v) == 1) v else NA }
cells <- do.call(rbind, lapply(names(CELLS), function(k) { cl <- CELLS[[k]]; x <- sets[sets$cell == k, ]
  g1 <- x[x$gen == "g1" & !is.na(x$recovered), ]; g2 <- x[x$gen == "g2" & !is.na(x$recovered), ]
  rd <- sapply(c(cl$g, cl$other), function(cv) realDIC(cl$task, cl$arch, cv))
  real21 <- rd[cl$g[2]] - rd[cl$g[1]]; real_b2 <- min(rd[-1]) - rd[cl$g[1]]
  s1 <- x$d21[x$gen == "g1" & !is.na(x$d21)]; s2 <- x$d21[x$gen == "g2" & !is.na(x$d21)]
  lr <- if (length(s1) > 1 && length(s2) > 1) dnorm(real21, mean(s2), sd(s2), log = TRUE) - dnorm(real21, mean(s1), sd(s1), log = TRUE) else NA
  data.frame(cell = k, task = cl$task, arch = cl$arch, G1 = cl$g[1], G2 = cl$g[2],
             sets_g1 = nrow(g1), recovered_g1 = sum(g1$recovered), false_alarms = sum(g1$winner_mech == 2),
             fa_max_margin = if (any(g1$winner_mech == 2)) max(g1$gap[g1$winner_mech == 2]) else 0,
             sets_g2 = nrow(g2), recovered_g2 = sum(g2$recovered), misses = sum(g2$winner == cl$g[1]),
             miss_max_margin = if (any(g2$winner == cl$g[1])) max(g2$gap[g2$winner == cl$g[1]]) else 0,
             floor = suppressWarnings(max(c(g1$gap, g2$gap))),
             d21_g1_mean = mean(s1), d21_g1_min = suppressWarnings(min(s1)), d21_g1_max = suppressWarnings(max(s1)),
             d21_g2_mean = mean(s2), d21_g2_min = suppressWarnings(min(s2)), d21_g2_max = suppressWarnings(max(s2)),
             real_d21 = real21, real_best2_minus_g1 = real_b2, log_lr_g2_vs_g1 = round(lr, 1),
             ll_flagged_sets = sum(x$ll_flag), row.names = NULL)
}))
write_out(cells, file.path(R5, "recovery_rec5_cells.csv"))
cat("fits with DIC:", sum(!is.na(map$DIC)), "of", nrow(map), "; stop rule met:", sum(map$met %in% TRUE),
    "; RDM fits screened:", sum(!is.na(map$ll_n_up)), "of", sum(map$arch == "RDM"), ", flagged:", sum(map$ll_n_up > 0, na.rm = TRUE), "\n\n")
print(cells[, c("cell", "task", "arch", "G1", "G2", "sets_g1", "recovered_g1", "false_alarms", "fa_max_margin", "sets_g2", "recovered_g2",
                "misses", "miss_max_margin", "floor")], row.names = FALSE, digits = 4)
cat("\nd21 = DIC(G2) - DIC(G1): simulated under each generator vs the real data (log_lr > 0 favours G2)\n")
print(cells[, c("cell", "d21_g1_mean", "d21_g1_min", "d21_g1_max", "d21_g2_mean", "d21_g2_min", "d21_g2_max", "real_d21", "log_lr_g2_vs_g1")],
      row.names = FALSE, digits = 4)
cat("\n"); print(sets[, c("cell", "gen", "replicate", "n_done", "n_met", "winner", "recovered", "gap", "d21", "ll_flag")], row.names = FALSE, digits = 4)
