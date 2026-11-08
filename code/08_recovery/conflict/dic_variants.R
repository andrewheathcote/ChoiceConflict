## Why every DIC in the project is the best-draw plug-in: recovery study 5 scored with both plug-ins.
##   Rscript code/08_recovery/conflict/dic_variants.R      (reads recovery/conflict/recovery_rec5.csv)
## Best-draw DIC: pD = meanD - min(minD, Dmean), DIC = meanD + pD. Standard DIC: pD = meanD - Dmean (deviance at the
## posterior mean), DIC = 2 meanD - Dmean. meanD, Dmean, minD are sums over subjects of the per-subject deviances; the
## best-draw minimum is taken per subject (EMC2 compare(use_best_fit = TRUE)), so the table's DIC is recomputed from the
## stored per-fit DIC rather than from the summed columns.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
m <- read.csv(file.path(ROOT, "recovery", "conflict", "recovery_rec5.csv"), stringsAsFactors = FALSE)
score <- function(v) do.call(rbind, lapply(split(m, m$set), function(x) { d <- x[[v]]; w <- which.min(d)
  data.frame(set = x$set[1], gen = x$gen[1], recovered = x$true_model[w], gap = d[x$true_model] - d[w]) }))
out <- do.call(rbind, lapply(c(best_draw = "DIC", standard = "DIC_standard"), function(v) { s <- score(v)
  data.frame(criterion = v, one_mech_recovered = sum(s$recovered[s$gen == "g1"]), one_mech_sets = sum(s$gen == "g1"),
             two_mech_recovered = sum(s$recovered[s$gen == "g2"]), two_mech_sets = sum(s$gen == "g2"),
             worst_wrong_winner_gap = round(max(s$gap))) }))
out$criterion <- c("best-draw plug-in (used)", "standard (posterior-mean plug-in)")
print(out, row.names = FALSE); write_out(out, file.path(ROOT, "recovery", "conflict", "dic_variants_rec5.csv"))
