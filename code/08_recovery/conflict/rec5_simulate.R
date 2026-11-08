## Recovery study 5, "one mechanism or two?": simulate the data sets and write the candidates' specs.
##   Rscript code/08_recovery/conflict/rec5_simulate.R      (needs the generators' fits in fits/, i.e. the local archive)
## Each cell pairs a one-mechanism generator (g1) with its two-mechanism rival (g2) on that cell's own task. The generator
## is the displayed fit (derived/tables/display_jm.csv); the data are simulated on the real task data's design: subject
## parameters drawn from MVN(posterior mean mu, posterior mean Sigma), then make_data(shrink2bound = TRUE) on the real
## subjects and trials. 7 cells x 2 generators x 5 replicates = 70 data sets; 4 candidates each = 280 fits. Each candidate
## is the real fit's own spec (design, prior, chain count) with only the data swapped.
## Writes, under $CC_RECOVERY_WORK (default <CC_ROOT>/recovery_work; not part of the archive):
##   inputs/data/<task>/r5_c<k>_g<j>_r<rep>.RData   simulated data        truth/truth_c<k>_g<j>_r<rep>.RData (alpha, mu, Sigma)
##   inputs/specs/r5_c<k>_g<j>_r<rep>_<arch>_<cand>.rds                     rec5_map.csv, cells.rds, sim_summary.csv
## Then fit each spec:  CC_INPUTS=<work>/inputs CC_REFITS=<work>/fits Rscript code/02_fit/fit.R <spec> jmh-rec5
## and score with score_rec5.R. From newsampler/rec5/build_rec5.R and make_specs_rec5.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R")); source(file.path(CODE, "lib", "predict.R"))
use_build(BUILD[["predict"]]); suppressMessages(library(mvtnorm))
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work")); REPS <- 5
for (d in c("inputs/specs", "truth")) dir.create(file.path(WORK, d), recursive = TRUE, showWarnings = FALSE)
CELLS <- list(
  c1 = list(task = "Simon",  arch = "RDM", g = c("Bi", "ViSi"), other = c("BiSi", "VcSi")),
  c2 = list(task = "Simon",  arch = "LBA", g = c("Bi", "Bci"),  other = c("BiSi", "ViSi")),
  c3 = list(task = "Stroop", arch = "LNR", g = c("Mi", "MiSi"), other = c("Mci", "MiSc")),
  c4 = list(task = "Stroop", arch = "LBA", g = c("Vi", "Vci"),  other = c("BcVi", "ViSi")),
  c5 = list(task = "ANT",    arch = "RDM", g = c("Bi", "BiVc"), other = c("BiVi", "ViSi")),
  c6 = list(task = "ANT",    arch = "LBA", g = c("Bi", "ViSi"), other = c("BiSi", "BiVi")),
  c7 = list(task = "Stroop", arch = "RDM", g = c("Vi", "ViSi"), other = c("Vci", "BiSi")))
saveRDS(CELLS, file.path(WORK, "cells.rds"))
D <- read_tab("display_jm.csv")
summ1 <- function(dat) { cong <- as.character(dat$RS) == as.character(dat$IS); ok <- as.character(dat$R) == as.character(dat$RS)
  qd <- function(k) mean(quantile(dat$rt[k & ok], c(.1, .9)))
  data.frame(acc_C = mean(ok[cong]), acc_I = mean(ok[!cong]), delta = qd(!cong) - qd(cong), mean_rt = mean(dat$rt)) }
S <- NULL; map <- NULL
for (k in names(CELLS)) { cl <- CELLS[[k]]
  real <- load_dat(cl$task); subs <- levels(factor(real$subjects))
  S <- rbind(S, cbind(cell = k, task = cl$task, generator = "real data", replicate = NA, summ1(real)))
  dir.create(file.path(WORK, "inputs/data", cl$task), recursive = TRUE, showWarnings = FALSE)
  for (j in 1:2) {
    spec <- paste(cl$task, cl$arch, cl$g[j], sep = "_"); lab <- D$label[D$spec == spec]; stopifnot(length(lab) == 1)
    f <- fit_path(lab); e <- load_fit(lab)
    des <- fit_design(f, real)
    mu <- rowMeans(get_pars(e, selection = "mu", return_mcmc = FALSE, merge_chains = TRUE))
    Sg <- apply(get_pars(e, selection = "Sigma", return_mcmc = FALSE, merge_chains = TRUE, remove_constants = FALSE), c(1, 2), mean)
    stopifnot(identical(names(mu), rownames(Sg))); rm(e); gc(verbose = FALSE)
    for (r in seq_len(REPS)) {
      set.seed(20261102 + 1000 * r + 10 * as.integer(sub("c", "", k)) + j)
      alpha <- rmvnorm(length(subs), mean = mu, sigma = Sg); rownames(alpha) <- subs; colnames(alpha) <- names(mu)
      dat <- make_data(alpha, design = des, data = real, shrink2bound = TRUE)
      stopifnot(nrow(dat) == nrow(real), identical(levels(factor(dat$subjects)), subs))
      tag <- sprintf("%s_g%d_r%d", k, j, r)
      save(dat, file = file.path(WORK, "inputs/data", cl$task, paste0("r5_", tag, ".RData")))
      save(alpha, mu, Sg, lab, file = file.path(WORK, "truth", paste0("truth_", tag, ".RData")))
      S <- rbind(S, cbind(cell = k, task = cl$task, generator = spec, replicate = r, summ1(dat)))
      for (cv in c(cl$g, cl$other)) {               # the four candidates' specs, data swapped
        sp <- read_spec(sprintf("%s_%s_%s", cl$task, cl$arch, cv)); nm <- sprintf("r5_%s_%s_%s", tag, cl$arch, cv)
        sp$data_file <- sprintf("r5_%s.RData", tag)
        sp$recovery <- list(study = "rec5", cell = k, task = cl$task, generator = spec, replicate = r,
                            candidate = paste(cl$task, cl$arch, cv, sep = "_"), true_model = cv == cl$g[j])
        saveRDS(sp, file.path(WORK, "inputs/specs", paste0(nm, ".rds")))
        map <- rbind(map, data.frame(spec = nm, cell = k, task = cl$task, arch = cl$arch, gen = paste0("g", j), generator = spec,
                                     replicate = r, candidate = cv, n_mech = if (cv == cl$g[1]) 1L else 2L,
                                     true_model = cv == cl$g[j], chains = sp$n_chains))
      }
    }
    cat(sprintf("%s g%d %-16s (%s): %d sets\n", k, j, spec, lab, REPS))
  }
}
stopifnot(!anyDuplicated(tolower(map$spec)), nrow(map) == 280)
write.csv(map, file.path(WORK, "rec5_map.csv"), row.names = FALSE); write.csv(S, file.path(WORK, "sim_summary.csv"), row.names = FALSE)
cat(nrow(map), "specs written to", file.path(WORK, "inputs/specs"), "\n")
