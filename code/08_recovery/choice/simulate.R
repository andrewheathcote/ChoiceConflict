## Choice recovery study (2026-11): is the DDM's DIC lead over the LBA on the simple-choice data recoverable?
##   Rscript code/08_recovery/choice/simulate.R     (needs fits/Choice/, i.e. the local archive)
## Generators: the choice fits Choice_DDM_v0_jmd1 and Choice_LBA_jmd1x (mean rate free), on their own build (4f1fa123). As the conflict
## study 5: subject parameters drawn from MVN(posterior mean mu, posterior mean Sigma), then make_data(shrink2bound =
## TRUE) on the real subjects and trials, design from the generator's spec. 5 replicates per generator. Candidates: the
## best fit of each family (DDM_v0, LBA, LNR, RDM, WDM_v0; the race models' mean rate free), each its own spec with only the data swapped:
## 2 x 5 x 5 = 50 fits. Writes, under $CC_RECOVERY_WORK (default <CC_ROOT>/recovery_work; not archived):
##   inputs/data/Choice/rc_<GEN>_r<r>.RData, truth/truth_rc_<GEN>_r<r>.RData, inputs/specs/rc_<GEN>_r<r>_<cand>.rds,
##   recch_map.csv, sim_summary_choice.csv
## Then fit each spec: CC_INPUTS=<work>/inputs CC_REFITS=<work>/fits Rscript code/02_fit/fit.R <spec> jmdrc2
## and score with dic.R, loo.R, score.R. From newsampler/gadi/choice_rec/sim_recch.R and make_specs_recch.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
use_build(BUILD[["choice"]]); suppressMessages(library(mvtnorm))
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work"))
for (d in c("inputs/specs", "inputs/data/Choice", "truth")) dir.create(file.path(WORK, d), recursive = TRUE, showWarnings = FALSE)
GEN <- c(DDM = "Choice_DDM_v0_jmd1", LBA = "Choice_LBA_jmd1x"); REPS <- 5      # generator fits (labels)
CAND <- c("DDM_v0", "LBA", "LNR", "RDM", "WDM_v0")
FT <- read_tab("FITS.csv"); gspec <- setNames(FT$spec[match(GEN, FT$label)], names(GEN)); stopifnot(!anyNA(gspec))
real <- load_dat("Choice"); subs <- levels(factor(real$subjects))
summ1 <- function(d) { ok <- as.character(d$R) == as.character(d$S)
  data.frame(acc = mean(ok), rt_mean = mean(d$rt), rt_q10 = quantile(d$rt[ok], .1), rt_q50 = quantile(d$rt[ok], .5),
             rt_q90 = quantile(d$rt[ok], .9), err_rt_mean = mean(d$rt[!ok]), row.names = NULL) }
S <- cbind(generator = "real data", replicate = NA, summ1(real)); map <- NULL
for (g in names(GEN)) {
  sp <- read_spec(gspec[[g]]); e <- load_fit(GEN[[g]])
  des <- design(data = real, model = get(sp$model), matchfun = sp$matchfun, functions = sp$funs, formula = sp$flist,
                contrasts = sp$clist, constants = sp$consts, pre_transform = sp$pre, report_p_vector = FALSE)
  mu <- rowMeans(get_pars(e, selection = "mu", return_mcmc = FALSE, merge_chains = TRUE))
  Sg <- apply(get_pars(e, selection = "Sigma", return_mcmc = FALSE, merge_chains = TRUE, remove_constants = FALSE), c(1, 2), mean)
  stopifnot(identical(names(mu), rownames(Sg)), identical(names(mu), sp$stored_par_names)); rm(e); gc(verbose = FALSE)
  for (r in seq_len(REPS)) {
    set.seed(20261104 + 1000 * r + match(g, names(GEN)))
    alpha <- rmvnorm(length(subs), mean = mu, sigma = Sg); rownames(alpha) <- subs; colnames(alpha) <- names(mu)
    dat <- make_data(alpha, design = des, data = real, shrink2bound = TRUE)
    stopifnot(nrow(dat) == nrow(real), identical(levels(factor(dat$subjects)), subs), !anyNA(dat$rt))
    tag <- sprintf("rc_%s_r%d", g, r)
    save(dat, file = file.path(WORK, "inputs/data/Choice", paste0(tag, ".RData")))
    save(alpha, mu, Sg, file = file.path(WORK, "truth", paste0("truth_", tag, ".RData")))
    S <- rbind(S, cbind(generator = gspec[[g]], replicate = r, summ1(dat)))
    for (cv in CAND) {
      cs <- read_spec(paste0("Choice_", cv)); nm <- sprintf("%s_%s", tag, cv)
      cs$data_file <- paste0(tag, ".RData")
      cs$recovery <- list(study = "choice recovery", generator = gspec[[g]], replicate = r, candidate = paste0("Choice_", cv),
                          true_model = cv == sub("^Choice_", "", gspec[[g]]))
      saveRDS(cs, file.path(WORK, "inputs/specs", paste0(nm, ".rds")))
      map <- rbind(map, data.frame(spec = nm, generator = gspec[[g]], replicate = r, candidate = paste0("Choice_", cv),
                                   family = sub("_.*", "", cv), true_model = cv == sub("^Choice_", "", gspec[[g]]), chains = cs$n_chains))
    }
  }
}
stopifnot(nrow(map) == 50, !anyDuplicated(tolower(map$spec)))
write.csv(map, file.path(WORK, "recch_map.csv"), row.names = FALSE); write.csv(S, file.path(WORK, "sim_summary_choice.csv"), row.names = FALSE)
print(S, digits = 3, row.names = FALSE)
