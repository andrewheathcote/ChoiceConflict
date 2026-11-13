## PSIS-LOO and WAIC for ONE fit, from one set of trialwise log-likelihoods, thinned to n_draws evenly over chains
## (DMC2 f2f58985: EMC2:::trialwise_pointwise, the code compare(LOO = TRUE) runs; r_eff from each subject's
## log-likelihood ESS). Used for the recovery fits (studies 1-5, choice recovery; on AWS) and the real fits
## (code/04_select/loo_real.R). Read-only on the fit.
##   Rscript loo.R <fit .RData> <out .rds> [n_draws = 2000] [ll_threads = 1]
## Build: EMC2_LIB (a library holding EMC2 f2f58985 and EMC2_SHA.txt); default <archive>/lib/f2f58985.
## Output <out .rds> (written via a temporary file, so its presence means done): summary (one row: elpd_loo, LOO =
## -2 elpd_loo, p_loo, se, mcse, k > 0.7 and > 1 among non-constant trials, WAIC and p_waic from the same pointwise
## values, draws used, minutes, seed, sha) and the pointwise vectors (subject, elpd_loo, pareto_k, constant, elpd_waic).
## The seed is fixed per fit (from its label); the computation itself draws no random numbers.
a <- commandArgs(TRUE); f <- a[1]; out <- a[2]
n_draws <- if (length(a) >= 3) as.integer(a[3]) else 2000L
threads <- if (length(a) >= 4) as.integer(a[4]) else 1L
LIB <- Sys.getenv("EMC2_LIB", "")
if (!nzchar(LIB)) {
  d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
  while (!file.exists(file.path(d, "config.R")) && d != dirname(d)) d <- dirname(d)
  LIB <- file.path(dirname(d), "lib", "f2f58985")
}
suppressMessages(library(EMC2, lib.loc = LIB)); SHA <- readLines(file.path(LIB, "EMC2_SHA.txt"))[1]
stopifnot(startsWith(SHA, "f2f58985"))
if (threads > 1) options(emc.ll_backend = "multithreaded", emc.n_threads = threads)
label <- sub("\\.RData$", "", basename(f))
seed <- 20261110L + sum(utf8ToInt(label)); set.seed(seed)
t0 <- Sys.time()
e <- new.env(); nm <- load(f, envir = e); emc <- e[[nm[1]]]; rm(e)
subs <- names(emc[[1]]$data)
pw <- EMC2:::trialwise_pointwise(emc, stage = "sample", filter = 0, subjects = subs, LOO = TRUE, r_eff = "subject",
                                 cores = 1, n_draws = n_draws)
loo <- pw$loo; waic <- pw$waic; const <- pw$constant
subj <- rep(subs, vapply(emc[[1]]$data[subs], function(d) EMC2:::n_trials_dadm(d), 1))
stopifnot(length(subj) == nrow(loo), length(const) == nrow(loo))
alpha1 <- get_pars(emc, selection = "alpha", stage = "sample", by_subject = TRUE, merge_chains = FALSE)[[1]]
draws <- sum(vapply(alpha1, function(ch) { n <- nrow(ch)
  length(unique(round(seq(1, n, length.out = min(n, ceiling(n_draws / length(alpha1))))))) }, 1))
k <- loo[, "influence_pareto_k"]; kk <- k[!const]
row <- data.frame(label = label, n_subjects = length(subs), n_trials = nrow(loo), n_constant = sum(const), draws = draws,
                  elpd_loo = sum(loo[, "elpd_loo"]), se_elpd_loo = sqrt(nrow(loo) * var(loo[, "elpd_loo"])),
                  LOO = -2 * sum(loo[, "elpd_loo"]), p_loo = sum(loo[, "p_loo"]),
                  mcse_elpd_loo = sqrt(sum(loo[, "mcse_elpd_loo"]^2)),
                  k_gt_0.7 = sum(kk > 0.7), k_gt_1 = sum(kk > 1), share_k_gt_0.7 = mean(kk > 0.7),
                  WAIC = -2 * sum(waic[, "elpd_waic"]), p_waic = sum(waic[, "p_waic"]),
                  minutes = as.numeric(difftime(Sys.time(), t0, units = "mins")), ll_threads = threads,
                  n_draws_arg = n_draws, seed = seed, emc2_sha = SHA)
print(row)
tmp <- paste0(out, ".tmp")
saveRDS(list(summary = row, subject = subj, elpd_loo = loo[, "elpd_loo"], pareto_k = k, constant = const,
             elpd_waic = waic[, "elpd_waic"], file = basename(f)), tmp)
invisible(file.rename(tmp, out))
