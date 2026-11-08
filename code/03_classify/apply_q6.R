## The display table with the Q6 refits (Andrew, 2026-11-07). The Simon LNR models with a congruency term on the
## location m (the Mc family: Mc, McSc, McSci, McSi, Mci, MciSc, MciSci, MciSi) were fitted with m_CIincong free,
## unlike ANT and Stroop; in Mci and MciSc it went the wrong way (incongruent faster). They were refit with the
## theory's order restriction m_CIincong = exp(.) (02_fit/fit_options.csv rows jmr and jmh1r):
##   jmr   : the jm settings on d0b6cf0b (Mci, MciSc, MciSi, Mc, McSc, McSci, McSi met the stop rule)
##   jmh1r : MciSci, which did not meet the rule on the jm settings (as its unrestricted jm fit), on the wave-2
##           campaign settings (9e8727dd, culling, flat test), as its unrestricted jmh1 fit was
## Each refit's class: converged when it met the stop rule (max(mu Rhat, q.99 of alpha Rhat) < 1.1) and no subject
## was flagged by the subject-mode screen (screen_q6.csv: T = 20 LL against models nested in the restricted space).
##   Rscript code/03_classify/apply_q6.R    display_jm_base.csv + summ_refit/ -> derived/tables/display_jm.csv
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
D  <- read_tab("display_jm_base.csv")
SC <- read_tab("screen_q6.csv")
Q6 <- c(Mc = "jmr", McSc = "jmr", McSci = "jmr", McSi = "jmr", Mci = "jmr", MciSc = "jmr", MciSci = "jmh1r", MciSi = "jmr")
for (v in names(Q6)) {
  spec <- paste0("Simon_LNR_", v); lab <- paste0(spec, "_", Q6[[v]])
  s <- readRDS(file.path(DERIVED, "summ_refit", paste0(lab, ".rds")))
  rule <- max(s$rhat[["mu"]], s$rhat[["alpha_q99"]])
  flagged <- SC$n_flagged[SC$stem == lab]; stopifnot(length(flagged) == 1)
  i <- which(D$spec == spec); stopifnot(length(i) == 1)
  D$label[i]  <- lab
  D$source[i] <- sprintf("%s (Q6 restricted refit)", Q6[[v]])
  D$rule[i]   <- rule
  D$DIC[i]    <- s$dic$DIC[1]
  D$class[i]  <- if (rule < 1.1 && flagged == 0) "converged" else if (rule < 1.1) "converged, subject-mode caveat" else "not converged"
  D$display_converged[i] <- D$class[i] == "converged"
  D$caveat[i] <- ""
}
write_out(D, tab_path("display_jm.csv"))
print(D[D$spec %in% paste0("Simon_LNR_", names(Q6)), c("spec", "label", "class", "rule", "DIC")], row.names = FALSE)
