## choice_jm.csv: one row per simple-choice fit (label jmd1), from its .txt record and choice_dic.csv
## (PLAN_choice_refit.md, "After the fits" 2). Kept separate from display_jm.csv.
##   Rscript code/04_select/choice_tables.R    (records in derived/records; writes derived/tables/choice_jm.csv)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
## class: "converged" when fit() met the stop rule (max_gd 1.1, gd_quantile .99, alpha and mu) in at least 1,000 sample
## draws, else "not converged". rule_all_draws = max(rhat_mu, rhat_alpha_q99) on all sample draws; fit()'s own check can
## drop the first third, so a converged fit can show a value just above 1.1 here.
txt <- list.files(file.path(DERIVED, "records"), "^Choice_.*_jmd1[.]txt$", full.names = TRUE)
rd <- function(f) { x <- readLines(f); k <- sub(":.*$", "", x); v <- sub("^[^:]*: ", "", x); v <- v[!duplicated(k, fromLast = TRUE)]
  setNames(v, k[!duplicated(k, fromLast = TRUE)]) }
ps <- function(s, re) as.numeric(sub(re, "\\1", regmatches(s, regexpr(re, s))))
R <- do.call(rbind, lapply(txt, function(f) { r <- rd(f); p <- r[["posthoc_stats"]]
  data.frame(spec = r[["spec"]], label = paste0(r[["spec"]], "_", r[["label"]]),
             family = sub("^(LBA|RDM|LNR|WDM|DDM).*$", "\\1", sub("^Choice_", "", r[["spec"]])),
             status = r[["status"]], stop_rule_met = as.logical(r[["stop_rule_met"]]),
             rhat_mu = as.numeric(r[["rhat_mu"]]), rhat_alpha_q99 = as.numeric(r[["rhat_alpha_q99"]]),
             rhat_alpha = as.numeric(r[["rhat_alpha"]]),
             sample_draws = as.integer(strsplit(r[["chain_n"]], "/")[[1]][4]), chains = as.integer(r[["chains"]]),
             ess_min = ps(p, "ess_min ([0-9.]+)"), ll_gap = ps(p, "ll_gap ([0-9.]+)"),
             subj_spread = ps(p, "subj_spread ([0-9.]+)"), subjects_spread_gt5 = ps(p, "spread > 5: ([0-9]+)"),
             posthoc_stats = p, minutes = as.numeric(r[["minutes"]]), max_tries = as.integer(r[["max_tries"]]),
             job = r[["job"]], emc2_sha = r[["emc2_sha"]], constraint = r[["constraint"]])
}))
R$rule_all_draws <- pmax(R$rhat_mu, R$rhat_alpha_q99)
R$class <- ifelse(R$status == "ok" & R$stop_rule_met, "converged", "not converged")
R$caveat <- ifelse(R$class == "converged", "",
                   sprintf("max_tries %d used; rule %.3f; accepted as is (Andrew 2026-10-23)", R$max_tries, R$rule_all_draws))
D <- read_tab("choice_dic.csv")
R <- merge(R, D[, c("spec", "variant", "n_par", "DIC", "dDIC", "Dmean", "meanD", "EffectiveN", "BPIC")], by = "spec")
R <- R[order(R$DIC), c("spec", "variant", "family", "label", "class", "caveat", "n_par", "DIC", "dDIC", "Dmean", "meanD",
                       "EffectiveN", "BPIC", "rule_all_draws", "rhat_mu", "rhat_alpha_q99", "rhat_alpha", "sample_draws",
                       "chains", "ess_min", "ll_gap", "subj_spread", "subjects_spread_gt5", "minutes", "max_tries", "job",
                       "emc2_sha", "constraint")]
write_out(R, tab_path("choice_jm.csv"))
print(R[, c("variant", "class", "n_par", "DIC", "dDIC", "Dmean", "rule_all_draws", "ess_min", "ll_gap", "subj_spread")],
      row.names = FALSE)
