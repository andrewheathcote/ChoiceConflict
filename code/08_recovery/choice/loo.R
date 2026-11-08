## WAIC, PSIS-LOO (pseudo-BMA+ with Bayesian bootstrap, stacking, k_bad) and best-draw DIC for one simulated data set of
## the choice recovery study: its five candidate fits, scored as 04_select/compare_choice.R scores the real fits (DMC2
## main f2f58985, use_best_fit = TRUE so DIC matches dic.R). Read-only on the fits.
##   Rscript code/08_recovery/choice/loo.R <generator: DDM|LBA> <replicate> [cores]   -> recovery/choice/loo/rc_<GEN>_r<r>.{csv,rds}
## From newsampler/gadi/choice_rec/recch_loo.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
use_build(BUILD[["ic"]])
a <- commandArgs(TRUE); G <- a[1]; R <- as.integer(a[2]); cores <- if (length(a) > 2) as.integer(a[3]) else 4L
set <- sprintf("rc_%s_r%d", G, R)
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work")); OUT <- file.path(ROOT, "recovery", "choice", "loo")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
map <- read.csv(file.path(WORK, "recch_map.csv"), stringsAsFactors = FALSE)
m <- map[startsWith(map$spec, paste0(set, "_")), ]; stopifnot(nrow(m) == 5)
sL <- setNames(lapply(m$spec, function(s) { e <- new.env(); nm <- load(file.path(WORK, "fits", paste0(s, "_jmdrc.RData")), envir = e); e[[nm[1]]] }), m$family)
set.seed(20261105 + R)
out <- compare(sL, BayesFactor = FALSE, WAIC = TRUE, LOO = TRUE, cores_loo = cores, use_best_fit = TRUE, print_summary = FALSE)
print(round(out, 3))
saveRDS(list(table = out, loo = attr(out, "loo"), set = set, emc2 = "f2f58985"), file.path(OUT, paste0(set, ".rds")))
write_out(cbind(set = set, generator = G, replicate = R, model = rownames(out), out), file.path(OUT, paste0(set, ".csv")))
