## Delta / CAF difference floors for the heat grids, for both group-curve aggregations (mean and median over subjects).
## Resample the subjects with replacement (2,000 times, seed 42), rebuild the observed group curve, and score it
## against the real observed curve with the SAME formula the models are scored with (06_predict/build_adequacy.R
## stat1: delta RMS over the five quantiles in % of the fixed effect size; CAF the mean of the incongruent and
## congruent RMS in % points); floor = 95th percentile. Group RT quantiles: mean or median over subjects; accuracy is
## always the mean. Observed primitives: lib/predict.R stats_one() on the task's data (Stroop: the accuracy-coded
## dat_acc.RData, as the original computation, which read them from a Stroop diffusion fit's summary).
## From newsampler/figs/noise_floor_agg.R.
##   Rscript code/05_floors/noise_floor_agg.R    -> derived/tables/noise_floors_agg.csv
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "predict.R"))
EFFECT <- c(Simon = 42, Stroop = 122, ANT = 70); B <- 2000
DATA <- c(Simon = "dat.RData", Stroop = "dat_acc.RData", ANT = "dat.RData")
med <- function(m) apply(m, 2, median, na.rm = TRUE); mn <- function(m) colMeans(m, na.rm = TRUE)
out <- NULL
for (tk in names(EFFECT)) {
  dat <- load_dat(tk, DATA[[tk]]); ob <- stats_one(dat, levels(factor(dat$subjects)))
  n <- nrow(ob$sqc)
  for (agg in c("mean", "median")) { af <- if (agg == "mean") mn else med
    O <- list(d = af(ob$sqi) - af(ob$sqc), aI = mn(ob$scafI), aC = mn(ob$scafC))
    set.seed(42); dl <- ca <- numeric(B)
    for (b in 1:B) { i <- sample.int(n, n, replace = TRUE)
      d <- af(ob$sqi[i, , drop = FALSE]) - af(ob$sqc[i, , drop = FALSE])
      dl[b] <- 100 * 1000 * sqrt(mean((d - O$d)^2)) / EFFECT[[tk]]
      ca[b] <- 100 * (sqrt(mean((mn(ob$scafI[i, , drop = FALSE]) - O$aI)^2)) + sqrt(mean((mn(ob$scafC[i, , drop = FALSE]) - O$aC)^2))) / 2 }
    out <- rbind(out, data.frame(measure = c("delta", "caf"), task = tk, agg = agg,
                                 floor = round(c(quantile(dl, .95), quantile(ca, .95)), 2), subjects = n)) }
}
write_out(out, tab_path("noise_floors_agg.csv")); print(out, row.names = FALSE)
