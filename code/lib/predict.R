## Subject-level posterior predictives (the only level used anywhere in this project).
## For each of n_post evenly spaced posterior draws, every real subject is simulated from its own alpha draw on the
## fit's own design (read back from the fit: read_fit(), then design()), with make_data(shrink2bound = TRUE).
## A draw that fails to simulate is left NA (the parameters are the subject's own, so it cannot be redrawn).
##   read_fit(f)              the design stored in a fit (formulas, contrasts, constants, functions, pre-transforms)
##   stats_one(d, subs)       per-subject correct-RT quantiles .1-.9 by congruency and five-bin accuracy (CAF)
##   pp_stats_subj(f, dat)    stats_one() for every draw: the primitives behind the heat grids' misfit
##   sim_raw(f, dat)          the simulated trials themselves (for EMC2's plot_delta / plot_caf)
## From NCI/rebuild_scripts.R (read_fit), NCI/adequacy_stats.R (stats_one) and newsampler/figs/pp_subj.R.
PROBS <- c(.1, .3, .5, .7, .9)
NBIN  <- 5

read_fit <- function(f) {
  ev <- new.env(); nm <- load(f, envir = ev); e <- ev[[nm[1]]]
  d  <- get_design(e)[[1]]
  list(variant = sub("[.]RData$", "", basename(f)), flist = d$Flist, clist = d$Clist, consts = d$constants,
       funs = d$Ffunctions, matchfun = d$matchfun, rlev = d$Rlevels, model = d$model()$c_name,
       pre = environment(d$model)$pre_transform, n_chains = length(e),
       npar = length(get_pars(e, selection = "mu", merge_chains = TRUE)$mu[[1]][1, ]))
}

## `subs` fixes the subject set and row order, so every draw returns matrices of the same shape even when a
## simulated subject has too few usable trials.
stats_one <- function(d, subs = NULL) {
  d <- d[is.finite(d$rt), ]
  cong <- d$RS == d$IS
  ## accuracy-coded data (Stroop WDM/DDM) carries R as FALSE/TRUE: correctness is the response itself
  corr <- if (identical(levels(factor(d$R)), c("FALSE", "TRUE"))) d$R == "TRUE" else d$R == d$RS
  if (is.null(subs)) subs <- levels(factor(d$subjects))
  subj <- factor(as.character(d$subjects), levels = subs)
  qs <- function(idx) {                        # per-subject correct-RT quantiles (at least 5 trials)
    m <- matrix(NA_real_, length(subs), length(PROBS), dimnames = list(subs, NULL))
    sp <- split(d$rt[idx], subj[idx])
    for (i in seq_along(sp)) if (length(sp[[i]]) >= 5) m[i, ] <- quantile(sp[[i]], PROBS, names = FALSE, type = 7)
    m
  }
  qc <- qs(cong & corr); qi <- qs(!cong & corr)
  caf <- function(idx) {                       # per-subject accuracy in NBIN equal-count RT bins (at least 2 x NBIN trials)
    m <- matrix(NA_real_, length(subs), NBIN, dimnames = list(subs, NULL))
    rs <- split(d$rt[idx], subj[idx]); cs <- split(corr[idx], subj[idx])
    for (i in seq_along(rs)) {
      r <- rs[[i]]; if (length(r) < NBIN * 2) next
      ## rank-based bins: near-degenerate predicted RTs make quantile breaks non-unique
      b <- ceiling(rank(r, ties.method = "first") / length(r) * NBIN); b[b > NBIN] <- NBIN
      m[i, ] <- tapply(cs[[i]], factor(b, levels = 1:NBIN), mean)
    }
    m
  }
  aI <- caf(!cong); aC <- caf(cong)
  list(gqc = colMeans(qc, na.rm = TRUE), gqi = colMeans(qi, na.rm = TRUE), sqc = qc, sqi = qi, scafI = aI, scafC = aC,
       gcafI = colMeans(aI, na.rm = TRUE), gcafC = colMeans(aC, na.rm = TRUE))
}

## the fit's design rebuilt on the data, and its subject-by-parameter draws
fit_design <- function(f, dat) {
  s <- read_fit(f)
  quietly(design(data = dat, model = get(s$model), matchfun = s$matchfun, functions = s$funs, formula = s$flist,
                 contrasts = s$clist, constants = s$consts, pre_transform = s$pre))
}
subject_draws <- function(f, subs) {
  ev <- new.env(); nm <- load(f, envir = ev); e <- ev[[nm[1]]]; rm(ev)
  al <- get_pars(e, selection = "alpha", merge_chains = TRUE, by_subject = TRUE, return_mcmc = FALSE)
  stopifnot(setequal(dimnames(al)[[2]], subs))
  al[, subs, , drop = FALSE]
}

pp_stats_subj <- function(f, dat, n_post = 100, seed = 42) {
  set.seed(seed)
  des <- fit_design(f, dat); subs <- levels(factor(dat$subjects)); ns <- length(subs)
  al <- subject_draws(f, subs); gc(verbose = FALSE)
  idx <- round(seq(1, dim(al)[3], length.out = n_post))
  GC <- GI <- matrix(NA_real_, n_post, length(PROBS)); CI <- CC <- matrix(NA_real_, n_post, NBIN)
  SC <- SI <- array(NA_real_, c(n_post, ns, length(PROBS))); AI <- AC <- array(NA_real_, c(n_post, ns, NBIN))
  bad <- 0L
  for (k in seq_along(idx)) {
    sim <- try(make_data(t(al[, , idx[k]]), design = des, data = dat, shrink2bound = TRUE), silent = TRUE)
    if (inherits(sim, "try-error")) { bad <- bad + 1L; next }
    st <- stats_one(sim, subs)
    GC[k, ] <- st$gqc; GI[k, ] <- st$gqi; CI[k, ] <- st$gcafI; CC[k, ] <- st$gcafC
    SC[k, , ] <- st$sqc; SI[k, , ] <- st$sqi; AI[k, , ] <- st$scafI; AC[k, , ] <- st$scafC
  }
  list(gqc = GC, gqi = GI, gcafI = CI, gcafC = CC, sqc = SC, sqi = SI, scafI = AI, scafC = AC,
       n_retry = 0, n_bad = bad, draws_used = idx)
}

## the same draws as pp_stats_subj() (same seed, same calls in the same order), kept as trials: stacked make_data()
## frames with a postn column, as EMC2's plot_delta / plot_caf take them. `fix` (named sampled values, e.g.
## c(v_iMd = -Inf)) overrides those parameters in every subject and draw, for predictives with a mechanism removed
## (06_predict/pp_null.R); same seed and draws, so only the fixed parameters differ from the fit's own predictives.
sim_raw <- function(f, dat, n_post = 100, seed = 42, fix = NULL) {
  set.seed(seed)
  des <- fit_design(f, dat); subs <- levels(factor(dat$subjects))
  al <- subject_draws(f, subs); gc(verbose = FALSE)
  stopifnot(all(names(fix) %in% dimnames(al)[[1]])); for (p in names(fix)) al[p, , ] <- fix[[p]]
  idx <- round(seq(1, dim(al)[3], length.out = n_post)); out <- vector("list", n_post)
  for (k in seq_along(idx)) {
    sim <- try(make_data(t(al[, , idx[k]]), design = des, data = dat, shrink2bound = TRUE), silent = TRUE)
    if (inherits(sim, "try-error")) next
    sim$postn <- k; out[[k]] <- sim
  }
  pp <- do.call(rbind, out[!vapply(out, is.null, TRUE)])
  pp <- pp[, c("postn", setdiff(names(pp), "postn"))]
  attr(pp, "n_retry") <- 0L; attr(pp, "n_draw") <- length(unique(pp$postn)); attr(pp, "draws_used") <- idx
  pp
}
## largest difference between stats_one() on each raw draw and stored primitives for that draw (a consistency check)
raw_vs_stored <- function(pp, stored, subs) {
  d <- 0
  for (k in unique(pp$postn)) { st <- stats_one(pp[pp$postn == k, ], subs)
    for (nm in c("sqc", "sqi", "scafI", "scafC")) d <- max(d, abs(st[[nm]] - stored[[nm]][k, , ]), na.rm = TRUE) }
  if (!identical(sort(unique(pp$postn)), which(!is.na(stored$gqc[, 1])))) d <- Inf
  d
}
