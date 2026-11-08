## Fit one model locally, with the settings of one of the project's fit sets (fit_options.csv), on this machine.
##   Rscript code/02_fit/fit.R <spec> <set> [cores] [chains]
##     spec   inputs/specs/<spec>.rds: design (formulas, contrasts, constants, functions, order restrictions), prior,
##            chain count and data file of the model
##     set    a row of fit_options.csv: the EMC2 build and the stop / culling / thinning settings that produced the
##            archived fits with that label suffix (e.g. jm, jmh-w2, jmd1, hf1)
##     cores  total cores (default: all detected); chains run in parallel, cores_per_chain = floor(cores / chains).
##            Race models use about one core per chain whatever cores_per_chain is; diffusion models scale within a chain.
## Writes <CC_REFITS, default refits/>/<spec>_<suffix>.RData and .txt (settings, timing, stop rule, culls); a second
## call with the same arguments resumes the partial file (sample checks counted over the whole fit); a finished fit is
## never overwritten. The seed is a function of the spec name and the set's seed label, as in every original run, so
## a refit starts where the original did; it reproduces the original fit statistically, not draw for draw (parallel
## chains, a different BLAS and the originals' resumes after interruptions all change the random path).
## Merged from newsampler/gadi/fit_new.R (jm), aws/jmh/fit_jmh.R (the campaign waves), gadi/choice/fit_dev.R (choice)
## and hardcases_niek/refit_local.R (hardfit); those originals are in local_only/historical/.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "transplant.R"))
a <- commandArgs(TRUE); if (length(a) < 2) stop("usage: fit.R <spec> <set> [cores] [chains]")
spec_name <- a[1]; set_name <- a[2]
cores <- if (length(a) > 2) as.integer(a[3]) else parallel::detectCores()
O <- read.csv(file.path(CODE, "02_fit", "fit_options.csv"), stringsAsFactors = FALSE, na.strings = c("", "NA"))
o <- O[O$set == set_name, ]; if (nrow(o) != 1) stop("no fit set ", set_name, " in fit_options.csv")
num <- function(k) if (is.na(o[[k]])) NULL else as.numeric(o[[k]])
SHA <- use_build(o$build)

sp <- read_spec(spec_name); dat <- spec_dat(sp)
is_ddm <- sp$arch %in% c("DDM", "WDM")
## chains: the 4th argument (FITS.csv records each archived fit's count: the campaign ran its DDMs with 3 or 6),
## else the set's, else the spec's
nc <- if (length(a) > 3) as.integer(a[4]) else if (is.na(o$chains)) sp$n_chains else as.integer(o$chains); sp$n_chains <- nc
cpc <- max(1, floor(cores / nc))
MAXT <- as.integer(if (is_ddm && !is.na(o$max_tries_ddm)) o$max_tries_ddm else o$max_tries)
seed_label <- if (is.na(o$seed_label)) o$suffix else o$seed_label
seed <- sum(utf8ToInt(spec_name)) * 100 + 10 + sum(utf8ToInt(seed_label)) %% 89

OUT <- Sys.getenv("CC_REFITS", file.path(ROOT, "refits")); dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(OUT, paste0(spec_name, "_", o$suffix)); txt <- paste0(stem, ".txt")
resume <- file.exists(paste0(stem, ".RData"))
if (resume && file.exists(txt) && any(grepl("^status: ok", readLines(txt)))) stop(stem, " is finished: never overwrite a fit")

## stop criteria: burn = EMC2's default rule (mean_gd 1.1 on alpha and mu), restated (fit() drops a stage's defaults
## once any criterion is given), with culling where the set had it; sample = the stop rule, max(mu Rhat, .99 quantile
## of the subject Rhats) < 1.1, plus the set's culling, flat test and opt-in criteria
crit <- list(burn = list(mean_gd = 1.1, omit_mpsrf = TRUE, selection = c("alpha", "mu")),
             sample = list(max_gd = 1.1, gd_quantile = .99, omit_mpsrf = TRUE, selection = c("alpha", "mu")))
if (isTRUE(o$burn_cull)) crit$burn$cull <- TRUE
if (!is.na(o$cull_chain_z_sample)) crit$sample$cull <- list(chain_z = o$cull_chain_z_sample)
if (!is.na(o$max_flat_loc)) crit$sample <- c(crit$sample, list(max_flat_loc = o$max_flat_loc, flat_mcse = num("flat_mcse")))
for (k in c("min_ess", "max_ll_gap", "max_subj_spread", "abort_gd", "abort_try")) if (!is.na(o[[k]])) crit$sample[[k]] <- o[[k]]

info <- c(spec = spec_name, set = set_name, label = paste0(spec_name, "_", o$suffix), seed_label = seed_label, seed = seed,
          chains = nc, cores_per_chain = cpc, total_cores = cores, rt_resolution = "1/60", emc2_sha = SHA,
          driver = o$driver, max_tries = MAXT, thin_target = o$thin_target,
          burn_criteria = paste(names(crit$burn), vapply(crit$burn, function(v) paste(unlist(v), collapse = "/"), ""), collapse = ", "),
          sample_criteria = paste(names(crit$sample), vapply(crit$sample, function(v) paste(unlist(v), collapse = "/"), ""), collapse = ", "),
          R = R.version.string, host = Sys.info()[["nodename"]], start = format(Sys.time()))
note <- function(k, v) { info[k] <<- as.character(v); writeLines(paste(names(info), info, sep = ": "), txt) }
note("start", info[["start"]])

args <- list(fileName = paste0(stem, ".RData"), cores_per_chain = cpc, cores_for_chains = nc)
## CC_SMOKE=k: build the model from the spec and run k preburn iterations, then stop (nothing is saved)
if (nzchar(Sys.getenv("CC_SMOKE"))) {
  set.seed(seed); emc <- build_emc(sp, dat, 1/60); t1 <- Sys.time()
  emc <- EMC2:::run_emc(emc, stage = "preburn", stop_criteria = list(iter = as.integer(Sys.getenv("CC_SMOKE"))),
                        cores_per_chain = cpc, cores_for_chains = nc)
  cat(sprintf("smoke: %s built (%d parameters, %d chains) and ran %s preburn iterations in %.1f s on %s\n", spec_name,
              length(emc[[1]]$par_names), length(emc), Sys.getenv("CC_SMOKE"), as.numeric(difftime(Sys.time(), t1, units = "secs")), SHA))
  unlink(txt); quit(save = "no")
}
t0 <- Sys.time()
res <- tryCatch({
  if (o$driver == "hardfit") {
    ## DMC2 main's hardfit(): the campaign stop rule, chain replacement on, early stop 1.5 at check 20, history proposals
    stopifnot(!resume, exists("hardfit"))
    set.seed(seed); emc <- build_emc(sp, dat, 1/60)
    emc <- init_chains(emc, particles = as.integer(o$init_particles), cores_for_chains = nc)
    r <- hardfit(emc, iter = 1000, step_size = 100, max_tries = MAXT, max_gd = 1.1, gd_quantile = .99, min_ess = NULL,
                 max_ll_gap = o$max_ll_gap, max_subj_spread = o$max_subj_spread, cull = TRUE, abort_gd = o$abort_gd,
                 abort_try = o$abort_try, cores_per_chain = 1, cores_for_chains = nc, fileName = args$fileName,
                 history_proposals = isTRUE(o$history_proposals), verbose = TRUE)
    emc <- r; save(emc, file = args$fileName); r
  } else {
    MT_RUN <- MAXT
    if (resume) {
      ev <- new.env(); nm <- load(paste0(stem, ".RData"), envir = ev); emc <- ev[[nm[1]]]; rm(ev)
      stopifnot(identical(emc[[1]]$par_names, sp$stored_par_names), length(emc) == nc)
      note("resume_from_chain_n", paste(chain_n(emc)[1, ], collapse = "/"))
      if (chain_n(emc)[1, "sample"] > 0) {       # max_tries (and abort_try) count sample checks over the whole fit
        done <- as.integer(chain_n(emc)[1, "sample"] %/% 100); MT_RUN <- max(0L, MAXT - done)
        if (!is.null(crit$sample$abort_try)) crit$sample$abort_try <- max(1L, crit$sample$abort_try - done)
        note("sample_checks_before_resume", done)
      }
    } else { set.seed(seed); emc <- build_emc(sp, dat, 1/60) }
    do.call(fit, c(list(emc, stop_criteria = crit, max_tries = MT_RUN), args))
  }
}, error = function(e) e)
ok <- !inherits(res, "error")
note("minutes", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2))
note("status", if (ok) "ok" else paste("error:", conditionMessage(res)))
if (ok) {
  note("chain_n", paste(chain_n(res)[1, ], collapse = "/"))
  g <- function(sel, q = 1) { v <- unlist(gd_summary(res, selection = sel, stage = "sample")); v <- v[is.finite(v)]
    if (q == 1) max(v) else unname(quantile(v, q)) }
  rmu <- g("mu"); raq <- g("alpha", .99)
  note("rhat_mu", round(rmu, 3)); note("rhat_alpha", round(g("alpha"), 3)); note("rhat_alpha_q99", round(raq, 3))
  note("stop_rule_value", round(max(rmu, raq), 3))
  note("stop_rule_met", max(rmu, raq) < 1.1 && chain_n(res)[1, "sample"] >= 1000)
  cu <- res[[1]]$culls; note("culls", length(cu))
  if (length(cu)) note("burn_exit_replacements", sum(vapply(cu, function(x) isTRUE(x$exit), TRUE)))
  if (!is.na(o$thin_target) && o$thin_target > 0 && exists("thin_to_ess")) {
    ## the campaign kept the unthinned object beside the thinned one
    full <- paste0(stem, "_full.RData"); stopifnot(file.copy(paste0(stem, ".RData"), full, overwrite = FALSE))
    emc <- thin_to_ess(res, target = o$thin_target); th <- emc[[1]]$thinning
    note("thinning", sprintf("interval %.3f, draws %d -> %d per chain", th$interval, as.integer(th$n_before), as.integer(th$n_after)))
    tmp <- paste0(stem, ".RData.tmp"); save(emc, file = tmp); stopifnot(file.rename(tmp, paste0(stem, ".RData")))
  }
}
note("end", format(Sys.time()))
cat(paste(names(info), info, sep = ": "), sep = "\n")
