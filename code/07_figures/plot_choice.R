## Figure 4 (simple choice) on the jmd1 refits: Paper.R's Choice block, redrawn from the new fits with
## subject-level posterior predictives (Andrew, 2026-10-23), with EMC2's plot_cdf / plot_delta / plot_caf as in the
## archive. Same layout as Paper.R's choice.pdf:
##   DIC bars for all 16 variants; defective CDFs for the five full variants (bias + unequal SDs); delta and CAF by
##   stimulus for the DIC-best full variant (BEST, default DDM; Paper.R drew the LBA, then the best)
##   Rscript code/07_figures/plot_choice.R sim [cores]    derived/pp_raw/Choice_<spec>_jmd1_subj.rds for the five (skips existing)
##   Rscript code/07_figures/plot_choice.R plot           figures/choice.pdf
## Predictives: for 100 evenly spaced posterior draws, every subject is simulated from its own alpha draw on the
## design rebuilt from the spec (as pp_subj.R raw, subj level). Simulated on the fits' own build (dev 4f1fa123).
## Inputs (read-only): fits/Choice/, inputs/specs/, inputs/data/Choice/, derived/tables/choice_dic.csv.
a <- commandArgs(TRUE); mode <- a[1]
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
OUT <- fig_path("choice.pdf")
FULL <- c(WDM = "Choice_WDM_v0", DDM = "Choice_DDM_v0", RDM = "Choice_RDM_v1", LNR = "Choice_LNR", LBA = "Choice_LBA_v1")
BEST <- Sys.getenv("BEST", "DDM")
raw_file <- function(spec) pp_raw_path(paste0(spec, "_jmd1"))
dat <- load_dat("Choice")

if (mode == "sim") {
  use_build(BUILD[["choice"]]); suppressMessages(library(parallel))
  cores <- if (length(a) > 1) as.integer(a[2]) else 8
  dir.create(file.path(DERIVED, "pp_raw"), showWarnings = FALSE)
  for (spec in FULL) {
    if (file.exists(raw_file(spec))) { cat(spec, ": exists\n"); next }
    t0 <- Sys.time(); sp <- read_spec(spec)
    des <- design(data = dat, model = get(sp$model), matchfun = sp$matchfun, functions = sp$funs, formula = sp$flist,
                  contrasts = sp$clist, constants = sp$consts, pre_transform = sp$pre, report_p_vector = FALSE)
    e <- load_fit(paste0(spec, "_jmd1"))
    stopifnot(identical(e[[1]]$par_names, sp$stored_par_names))
    al <- get_pars(e, selection = "alpha", stage = "sample", merge_chains = TRUE, by_subject = TRUE, return_mcmc = FALSE)
    subs <- levels(factor(dat$subjects)); stopifnot(setequal(dimnames(al)[[2]], subs)); al <- al[, subs, , drop = FALSE]
    rm(e); gc(verbose = FALSE)
    idx <- round(seq(1, dim(al)[3], length.out = 100))
    out <- mclapply(seq_along(idx), function(k) {
      set.seed(42 + k)
      sim <- try(make_data(t(al[, , idx[k]]), design = des, data = dat), silent = TRUE)
      if (inherits(sim, "try-error")) return(NULL)
      sim$postn <- k; sim }, mc.cores = cores)
    pp <- do.call(rbind, out[!vapply(out, is.null, TRUE)])
    attr(pp, "n_draw") <- length(unique(pp$postn)); attr(pp, "draws_used") <- idx
    saveRDS(pp, raw_file(spec))
    cat(sprintf("%s: %d of 100 draws, %d rows, %.1f min\n", spec, attr(pp, "n_draw"), nrow(pp),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
}

if (mode == "plot") {
  use_build(BUILD[["choice"]])
  source(file.path(CODE, "lib", "plot_compare.R"))            # plot_compare(), from the legacy Functions.R (as Paper.R)
  D <- read_tab("choice_dic.csv")
  ict <- data.frame(DIC = D$DIC, BPIC = D$BPIC, row.names = D$variant)   # two columns: plot_compare subsets rows
  pp <- lapply(FULL, function(s) readRDS(raw_file(s)))
  pct <- c(10, 50, 90)
  cfun <- \(d)factor(d$S!=d$R,labels=c("Correct","Error"))
  csfun <- \(d)factor(paste0(d$S!=d$R,d$S))
  ttl <- c(WDM = "Wiener Diffusion Model (WDM)", DDM = "Diffusion Decision Model (DDM)", RDM = "Racing Diffusion Model (RDM)",
           LNR = "Lognormal Race (LNR)", LBA = "Linear Ballistic Accumulator (LBA)")
  pdf(OUT, width = 14, height = 7)
  par(mfrow = c(2, 4))
  plot_compare(ict, ylim = c(0, 8000), do = c("WDM", "DDM", "RDM", "LBA", "LNR"), spaceby = 1)
  abline(h = ict["DDM", "DIC"] - min(ict$DIC), lty = 3)
  for (m in names(FULL))
    plot_cdf(dat, pp[[m]], defective_factor = "AccCon",
             functions = list(Accuracy = cfun, AccCon = csfun), layout = NULL, legendpos = c("topleft", NA), xlab = "RT (sec)",
             add_percentiles = pct, ylab = "Defective CDF",
             line_types = c(1, 2, 1, 2), panel_factor = "S", main = ttl[[m]],
             defective_legend = list(legend = c("1", "2", "Data", "Fit"),
                                     lty = c(1:2, NA, NA), pch = c(NA, NA, 16, 1), title = "Stimulus"))
  plot_delta(dat, pp[[BEST]], delta_factor = "Stimulus", main = paste(BEST, "Delta Function"), layout = NULL, legendpos = NA,
             add_percentiles = pct, functions = list(Stimulus = \(d)d$S), ylab = "RT (2-1) (sec)")
  plot_caf(dat, pp[[BEST]], caf_factor = "Stimulus", main = paste(BEST, "Conditional Accuracy Function"), layout = NULL,
           ylab = "Accuracy (%)", smooth_window = 10, legendpos = c("bottomright", NA), functions = list(Stimulus = \(d)d$S),
           accuracy_function = function(d) d$S == d$R, ylim = c(90, 100))
  invisible(dev.off())
  cat("wrote", OUT, "\n")
}
