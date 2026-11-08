## Reduce one fit to the quantities the classification and model-selection tables use. Read-only on the fit.
##   Rscript code/02_fit/summarise_fit.R <label> [build]     -> derived/summ/<label>.rds
##     label  an archived fit (fits/<task>/<arch>/<label>.RData) or a refit (refits/<label>.RData)
##     build  the EMC2 build to load (default: the fit's own, from derived/tables/FITS.csv)
## Contents:
##   q     per quantity on the sampled scale (mu, lsd = log sqrt diag Sigma, alpha): mean and SD of the sample-stage
##         draws, all chains pooled
##   rhat  largest split-Rhat of mu, lsd and alpha, and the .99 quantile of the alpha Rhats; the stop rule is
##         max(rhat["mu"], rhat["alpha_q99"]) < 1.1
##   ess   coda ESS summed over chains: smallest for mu and lsd, median and 5% quantile for alpha
##   dic   compare(..., BayesFactor = FALSE, WAIC = FALSE, type = "conditional", use_best_fit = TRUE): the conditional
##         DIC with the best-draw plug-in, pD = sum_s [meanD_s - min(minD_s, Dmean_s)] (every DIC in the project)
## Predictives are not made here: they are subject level only (06_predict/pp_subj.R).
## From newsampler/gadi/summarise_fit.R (08354987), which also stored hyper-level predictive primitives.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
a <- commandArgs(TRUE); label <- a[1]
FT <- read_tab("FITS.csv")
build <- if (length(a) > 1) a[2] else FT$build[FT$label == label]
if (length(build) != 1 || is.na(build)) stop(label, " is not in FITS.csv: give its build as the second argument")
if (grepl("^Niek", build)) build <- "61e21587"     # Niek's fits were scored on EMC2 dev 61e21587 (the RDM likelihood fixed)
SHA <- use_build(build)
f <- tryCatch(fit_path(label), error = function(e) file.path(Sys.getenv("CC_REFITS", file.path(ROOT, "refits")), paste0(label, ".RData")))
ev <- new.env(); nm <- load(f, envir = ev); e <- ev[[nm[1]]]; rm(ev); class(e) <- "emc"
t0 <- Sys.time()
s <- lapply(e, function(ch) ch$samples); keep <- s[[1]]$stage == "sample"
mu  <- simplify2array(lapply(s, function(x) t(x$theta_mu[, keep, drop = FALSE])))            # iter x par x chain
lsd <- simplify2array(lapply(s, function(x) t(apply(x$theta_var[, , keep, drop = FALSE], 3, function(V) 0.5 * log(diag(V))))))
al  <- simplify2array(lapply(s, function(x) { A <- x$alpha[, , keep, drop = FALSE]
  matrix(aperm(A, c(3, 1, 2)), nrow = sum(keep), dimnames = list(NULL, as.vector(outer(dimnames(A)[[1]], dimnames(A)[[2]], paste, sep = "|")))) }))
dimnames(lsd)[[2]] <- dimnames(mu)[[2]] <- rownames(s[[1]]$theta_mu)
pool <- function(X) apply(X, 2, c)
qrow <- function(cl, X) { P <- pool(X); data.frame(class = cl, name = colnames(P), mean = colMeans(P), sd = apply(P, 2, sd)) }
ess  <- function(X) apply(X, 2, function(v) sum(apply(v, 2, function(ch) coda::effectiveSize(ch))))
ra <- EMC2:::split_rhat(al)
out <- list(label = label, file = basename(f), emc2_sha = SHA, q = rbind(qrow("mu", mu), qrow("lsd", lsd), qrow("alpha", al)),
            draws = sum(keep), chains = length(e), chain_n = chain_n(e),
            rhat = c(mu = max(EMC2:::split_rhat(mu), na.rm = TRUE), lsd = max(EMC2:::split_rhat(lsd), na.rm = TRUE),
                     alpha = max(ra, na.rm = TRUE), alpha_q99 = unname(quantile(ra, .99, na.rm = TRUE))),
            ess = c(mu_min = min(ess(mu)), lsd_min = min(ess(lsd)), alpha_median = median(ess(al)), alpha_q05 = unname(quantile(ess(al), .05))))
ic <- try(compare(list(m = e), BayesFactor = FALSE, WAIC = FALSE, type = "conditional", use_best_fit = TRUE, print_summary = FALSE), silent = TRUE)
if (inherits(ic, "try-error"))                     # builds whose compare() has no use_best_fit argument: it is the default there
  ic <- try(compare(list(m = e), BayesFactor = FALSE, WAIC = FALSE, type = "conditional", print_summary = FALSE), silent = TRUE)
out$dic <- if (inherits(ic, "try-error")) as.character(ic) else ic
out$summary_minutes <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
dir.create(file.path(DERIVED, "summ"), recursive = TRUE, showWarnings = FALSE)
saveRDS(out, file.path(DERIVED, "summ", paste0(label, ".rds")))
cat(sprintf("%s: draws %d, rule %.3f, DIC %s, %.1f min\n", label, out$draws, max(out$rhat[["mu"]], out$rhat[["alpha_q99"]]),
            if (is.data.frame(out$dic)) format(round(out$dic$DIC[1], 1)) else "ERR", out$summary_minutes))
