## Subject-level posterior predictives for archived fits (lib/predict.R): 100 evenly spaced draws, seed 42, every
## real subject simulated from its own alpha draw on the fit's own design. On the predictive build (BUILD["predict"],
## 3c1bcaf1), which reproduced the campaign's stored predictive primitives to 2e-14 on fits from two other builds.
##   Rscript code/06_predict/pp_subj.R run <cores> [label ...]   derived/pp_subj/<label>.rds  (default: every displayed-
##                                                               converged model without one; existing files are kept)
##   Rscript code/06_predict/pp_subj.R raw <cores> <label> ...    derived/pp_raw/<label>_subj.rds: the simulated trials
##                                                               (for EMC2's plot_delta / plot_caf), checked against the
##                                                               primitives in pp_subj/<label>.rds
## Memory: up to ~11 GB per process for the largest ANT fits; cores = processes (one fit each).
## From newsampler/figs/pp_subj.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "predict.R"))
SHA <- use_build(BUILD[["predict"]]); suppressMessages(library(parallel))
FT <- read_tab("FITS.csv")
## EMC2's suppress_output() (R/model_DDM.R, also in DMC2 main) opens a tempfile() per unique parameter row of a WDM/DDM
## simulation and never removes it: ~1e5 files per fit, and R's exit clean-up of millions of them takes hours.
## Removed after each fit (results unaffected).
drop_sink_files <- function() unlink(list.files(tempdir(), "^file", full.names = TRUE))
spec_of <- function(label) FT$spec[FT$label == label]
a <- commandArgs(TRUE); mode <- a[1]; cores <- as.integer(a[2]); labels <- a[-(1:2)]
if (mode == "run") {
  out_dir <- file.path(DERIVED, "pp_subj"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (!length(labels)) { D <- read_tab("display_jm.csv"); labels <- D$label[D$display_converged] }
  todo <- labels[!file.exists(pp_subj_path(labels))]
  cat(length(todo), "of", length(labels), "to do on", cores, "cores, EMC2", substr(SHA, 1, 8), "\n")
  invisible(mclapply(todo, function(lb) {
    t0 <- proc.time()[["elapsed"]]; sp <- read_spec(spec_of(lb)); dat <- spec_dat(sp)
    pp <- try(pp_stats_subj(fit_path(lb), dat), silent = TRUE)
    if (inherits(pp, "try-error")) { cat(sprintf("  FAILED %-26s %s\n", lb, substr(gsub("[\r\n]+", " ", as.character(pp)), 1, 120))); return(NULL) }
    out <- list(label = lb, spec = spec_of(lb), level = "subj", pp = pp, obs = stats_one(dat, levels(factor(dat$subjects))),
                emc2_sha = SHA, seed = 42, minutes = (proc.time()[["elapsed"]] - t0) / 60)
    tmp <- file.path(out_dir, paste0(".", lb, ".rds")); saveRDS(out, tmp); file.rename(tmp, pp_subj_path(lb))
    cat(sprintf("  %-26s %5.1f min  bad %d\n", lb, out$minutes, pp$n_bad)); drop_sink_files(); gc(verbose = FALSE); NULL
  }, mc.cores = cores, mc.preschedule = FALSE))
} else if (mode == "raw") {
  out_dir <- file.path(DERIVED, "pp_raw"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  todo <- labels[!file.exists(pp_raw_path(labels))]
  cat(length(todo), "of", length(labels), "raw to do on", cores, "cores, EMC2", substr(SHA, 1, 8), "\n")
  invisible(mclapply(todo, function(lb) {
    t0 <- proc.time()[["elapsed"]]; sp <- read_spec(spec_of(lb)); dat <- spec_dat(sp)
    pp <- try(sim_raw(fit_path(lb), dat), silent = TRUE)
    if (inherits(pp, "try-error")) { cat(sprintf("  FAILED %-26s %s\n", lb, substr(gsub("[\r\n]+", " ", as.character(pp)), 1, 120))); return(NULL) }
    d <- raw_vs_stored(pp, readRDS(pp_subj_path(lb))$pp, levels(factor(dat$subjects)))
    attr(pp, "label") <- lb; attr(pp, "emc2_sha") <- SHA; attr(pp, "seed") <- 42; attr(pp, "max_diff_stored") <- d
    tmp <- file.path(out_dir, paste0(".", lb, "_subj.rds")); saveRDS(pp, tmp); file.rename(tmp, pp_raw_path(lb))
    cat(sprintf("  %-26s %5.1f min  %8d rows  %3d draws  max |raw - stored| %.2g  %4.0f MB\n", lb,
                (proc.time()[["elapsed"]] - t0) / 60, nrow(pp), attr(pp, "n_draw"), d, file.size(pp_raw_path(lb)) / 1e6))
    rm(pp); drop_sink_files(); gc(verbose = FALSE); NULL
  }, mc.cores = cores, mc.preschedule = FALSE))
} else stop("usage: pp_subj.R run <cores> [label ...] | raw <cores> <label> ...")
