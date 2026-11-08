## The DIC difference floor from the model-recovery studies: the largest margin by which a wrong model beat the
## generating one (gap of the winner over the generator), over conflict recovery studies 1-4 (recovery_jm_sets.csv:
## LNR 149.4, LBA 86.8, RDM 81.6, cross-architecture 0) and study 5 (recovery_rec5_sets.csv: 52.8), rounded UP to
## the next 10: 150 (Andrew, 2026-11-07). The same 150 is the DROP threshold of the heat grids' "=" rule.
##   Rscript code/05_floors/recovery_floor.R    -> inputs/noise_floors.csv (dic rows; delta / CAF floors are in
##                                                 derived/tables/noise_floors_agg.csv, from noise_floor_agg.R)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
S1 <- read.csv(file.path(ROOT, "recovery", "conflict", "recovery_jm_sets.csv"), stringsAsFactors = FALSE)
S5 <- read.csv(file.path(ROOT, "recovery", "conflict", "recovery_rec5_sets.csv"), stringsAsFactors = FALSE)
g <- c(tapply(S1$gap, S1$study, max, na.rm = TRUE), rec5 = max(S5$gap, na.rm = TRUE))
raw <- max(g); fl <- ceiling(raw / 10) * 10
print(round(g, 1)); cat(sprintf("largest recovery gap %.1f -> DIC floor %d\n", raw, fl))
write_out(data.frame(measure = "dic", task = TASKS, floor = fl, largest_gap = round(raw, 2)), file.path(INPUTS, "noise_floors.csv"))
