## Full compare() table (DMC2 main f2f58985) for the best simple-choice fit of each family (choice_dic.csv):
## DIC/BPIC (conditional, posterior-mean deviance), marginal deviance by bridge sampling (MD, wMD), WAIC, PSIS-LOO
## with pseudo-BMA+ (Bayesian bootstrap) and stacking weights, k_bad. Read-only on the fits.
##   Rscript code/04_select/compare_choice.R      (writes derived/tables/compare_choice.{csv,rds})
## DIC: best-draw plug-in (use_best_fit = TRUE, as every DIC in this project; f2f58985 defaults to standard DIC).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
use_build(BUILD[["ic"]])
LAB <- c(DDM = "Choice_DDM_v0_jmd1", LBA = "Choice_LBA_v1_jmd1", LNR = "Choice_LNR_jmd1",
         RDM = "Choice_RDM_v1_jmd1", WDM = "Choice_WDM_v0_jmd1")
sL <- lapply(LAB, load_fit)
## Bridge sampling per model (compare(BayesFactor = TRUE) stops on the first model whose bridge sampler does not
## converge, 4 Nov): MD = -2 log marginal likelihood, NA where it did not converge; wMD over the models that did.
t0 <- Sys.time(); set.seed(20261104)
MD <- sapply(names(sL), function(m) { t1 <- Sys.time()
  r <- tryCatch(-2 * EMC2:::run_bridge_sampling(sL[[m]], stage = "sample", filter = 0, both_splits = FALSE,
                                                cores_for_props = 4, cores_per_prop = 1),
                error = function(e) { cat(m, "bridge sampling failed:", conditionMessage(e), "\n"); NA_real_ })
  cat(m, "MD", r, "(", round(as.numeric(difftime(Sys.time(), t1, units = "mins")), 1), "min )\n"); r })
out <- compare(sL, BayesFactor = FALSE, WAIC = TRUE, LOO = TRUE, cores_loo = 2, print_summary = FALSE, use_best_fit = TRUE)
ok <- !is.na(MD); wMD <- rep(NA_real_, length(MD)); wMD[ok] <- exp(-(MD[ok] - min(MD[ok])) / 2) / sum(exp(-(MD[ok] - min(MD[ok])) / 2))
out <- cbind(MD = MD[rownames(out)], wMD = setNames(wMD, names(MD))[rownames(out)], out)
print(round(out, 3))
cat("minutes:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "\n")
print(attr(out, "loo"))
saveRDS(list(table = out, loo = attr(out, "loo"), labels = LAB, emc2 = "f2f58985", use_best_fit = TRUE), tab_path("compare_choice.rds"))
write_out(cbind(model = rownames(out), label = LAB[rownames(out)], out), tab_path("compare_choice.csv"))
