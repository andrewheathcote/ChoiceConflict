## Put a stored fit's draws into a freshly built emc object (new build's design and model, chosen rt_resolution).
## Used (a) to score a stored legacy posterior on 1/60 data, so its DIC is comparable with 1/60 refits, and
## (b) to warm-start from a stored fit's last draws inside a 1/60 object. In memory only; nothing is saved here.
##   build_emc(sp, dat, rt_resolution)  rebuild the stored fit's model from its spec (checks names and prior)
##   transplant(fresh, leg)             fresh object carrying leg's samples, chain by chain
build_emc <- function(sp, dat, rt_resolution) {
  des <- design(data = dat, model = get(sp$model), matchfun = sp$matchfun, functions = sp$funs,
                formula = sp$flist, contrasts = sp$clist, constants = sp$consts, pre_transform = sp$pre,
                report_p_vector = FALSE)
  emc <- make_emc(dat, des, type = sp$type, n_chains = sp$n_chains, rt_resolution = rt_resolution)
  stopifnot(identical(emc[[1]]$par_names, sp$stored_par_names),
            identical(as.character(emc[[1]]$subjects), as.character(sp$stored_subjects)))
  pr <- emc[[1]]$prior; ps <- sp$prior
  if (is.null(ps)) return(emc)                     # composed cell (no stored fit): default prior, names checked above
  for (k in c("theta_mu_mean", "theta_mu_var", "v", "A"))
    if (!isTRUE(all.equal(unname(unclass(pr[[k]])), unname(unclass(ps[[k]])))))
      stop("prior differs from the stored fit in ", k)
  emc
}
transplant <- function(fresh, leg) {
  stopifnot(length(fresh) == length(leg), identical(fresh[[1]]$par_names, leg[[1]]$par_names),
            identical(as.character(fresh[[1]]$subjects), as.character(leg[[1]]$subjects)))
  for (i in seq_along(fresh)) fresh[[i]]$samples <- leg[[i]]$samples
  class(fresh) <- "emc"
  fresh
}

## IC() takes the mean deviance from the log-likelihoods stored at sampling time (samples$subj_ll) and only
## recomputes the plug-in deviance on the object's data. To score a posterior on other data, recompute
## subj_ll for every draw on that data with EMC2's own likelihood. Returns the object with subj_ll replaced.
rescore <- function(e) {
  for (i in seq_along(e)) {
    A <- e[[i]]$samples$alpha                       # par x subject x iteration
    ll <- e[[i]]$samples$subj_ll                    # subject x iteration
    for (sub in dimnames(A)[[2]])
      ll[sub, ] <- EMC2:::calc_ll_manager(t(A[, sub, ]), dadm = e[[i]]$data[[sub]], e[[i]]$model)
    e[[i]]$samples$subj_ll <- ll
  }
  e
}
