## Subject-level predictives with one mechanism removed (Andrew, 2026-11-13): sim_raw() (lib/predict.R) on a fit with
## the mechanism's sampled parameter fixed at its null for every subject and draw, so its natural-scale effect is 0.
## Same 100 draws and seed as the fit's own pp_raw, so the two differ only through the removed mechanism.
## Nulls: an order-restricted (exp) effect at -Inf (exp(-Inf) = 0), any other effect at 0 (s_iMd is on the log scale).
##   Rscript code/06_predict/pp_null.R <cores> <label> <par> [<label> <par> ...]
## writes derived/pp_null/<label>_<par>0.rds (existing files are kept).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "predict.R"))
SHA <- use_build(BUILD[["predict"]]); suppressMessages(library(parallel))
FT <- read_tab("FITS.csv")
pp_null_path <- function(label, par) file.path(DERIVED, "pp_null", sprintf("%s_%s0.rds", label, par))
a <- commandArgs(TRUE); cores <- as.integer(a[1]); jobs <- matrix(a[-1], ncol = 2, byrow = TRUE)
out_dir <- file.path(DERIVED, "pp_null"); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
jobs <- jobs[!file.exists(pp_null_path(jobs[, 1], jobs[, 2])), , drop = FALSE]
cat(nrow(jobs), "to do on", cores, "cores, EMC2", substr(SHA, 1, 8), "\n")
invisible(mclapply(seq_len(nrow(jobs)), function(k) {
  lb <- jobs[k, 1]; p <- jobs[k, 2]; t0 <- proc.time()[["elapsed"]]
  sp <- read_spec(FT$spec[FT$label == lb]); dat <- spec_dat(sp)
  null <- if (p %in% names(sp$pre$func) && sp$pre$func[[p]] == "exp") -Inf else 0
  pp <- try(sim_raw(fit_path(lb), dat, fix = setNames(null, p)), silent = TRUE)
  if (inherits(pp, "try-error")) { cat(sprintf("  FAILED %-22s %-6s %s\n", lb, p, substr(gsub("[\r\n]+", " ", as.character(pp)), 1, 120))); return(NULL) }
  attr(pp, "label") <- lb; attr(pp, "fixed") <- setNames(null, p); attr(pp, "emc2_sha") <- SHA; attr(pp, "seed") <- 42
  dst <- pp_null_path(lb, p); tmp <- file.path(out_dir, paste0(".", basename(dst))); saveRDS(pp, tmp); file.rename(tmp, dst)
  cat(sprintf("  %-22s %-6s = %4s  %5.1f min  %8d rows  %3d draws\n", lb, p, null,
              (proc.time()[["elapsed"]] - t0) / 60, nrow(pp), attr(pp, "n_draw")))
  rm(pp); gc(verbose = FALSE); NULL
}, mc.cores = cores, mc.preschedule = FALSE))
