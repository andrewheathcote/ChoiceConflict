## Score the choice recovery study by best-draw DIC (recch_dic.csv from dic.R): per simulated data set, the DIC winner
## among the five candidates, the true model's gap to the winner, and the DDM - LBA gap (LBA DIC minus DDM DIC; positive
## = DDM ahead), set against the real data's gap (derived/tables/choice_dic.csv). score_recch_ic.R adds WAIC and LOO.
##   Rscript code/08_recovery/choice/score.R     -> recovery/choice/recch_score.csv
## From newsampler/gadi/choice_rec/score_recch.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
RC <- file.path(ROOT, "recovery", "choice")
D <- read.csv(file.path(RC, "recch_dic.csv"), stringsAsFactors = FALSE)
real <- read_tab("choice_dic.csv")
real_gap <- real$DIC[real$spec == "Choice_LBA_v1"] - real$DIC[real$spec == "Choice_DDM_v0"]
D$set <- paste(D$generator, D$replicate)
S <- do.call(rbind, lapply(split(D, D$set), function(d) {
  d <- d[order(d$DIC), ]
  data.frame(generator = sub("^Choice_", "", d$generator[1]), r = d$replicate[1], winner = d$family[1],
             runner_up = d$family[2], margin = round(d$DIC[2] - d$DIC[1], 1),
             true_gap = round(d$DIC[d$true_model] - d$DIC[1], 1),
             DDM_minus_LBA = round(d$DIC[d$family == "LBA"] - d$DIC[d$family == "DDM"], 1),
             order = paste(d$family, collapse = " < "),
             unmet = paste(d$family[d$stop_rule_met != "TRUE"], collapse = ","))
}))
S <- S[order(S$generator, S$r), ]; rownames(S) <- NULL
print(S); cat(sprintf("\nreal data: LBA - DDM = %.1f (DDM ahead)\n", real_gap))
for (g in unique(S$generator)) { s <- S[S$generator == g, ]
  cat(sprintf("%s generated: recovered %d/%d; DDM - LBA gap range %.0f to %.0f\n", g, sum(s$winner == sub("_v[01]$", "", g)), nrow(s),
              min(s$DDM_minus_LBA), max(s$DDM_minus_LBA))) }
write_out(S, file.path(RC, "recch_score.csv"))
