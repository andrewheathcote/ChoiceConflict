## Score conflict recovery studies 1-4 (specs from studies1to4_specs.R, fits with set jmh-rec1).
##   Rscript code/08_recovery/conflict/score_studies1to4.R
##   reads <work>/rec_map.csv, <work>/summ/<spec>_jmhrec1.rds (02_fit/summarise_fit.R with CC_DERIVED=<work>) and the fit
##   records <work>/fits/<spec>_jmhrec1.txt. rXA_g2_r2_RDM_VciSci was refit as jmhrec1b (its jmhrec1 save and summary failed
##   on a full disk); a _jmhrec1b summary, where present, is used in place of _jmhrec1.
## Per data set (study x generator x replicate) with every candidate fitted: the conditional-DIC winner, whether it is the
## generator, and the true-vs-winner gap, DIC(true) - DIC(winner) (0 when recovered). meanD is kept for bpic.R. The DIC floor is the largest gap
## (05_floors/recovery_floor.R). Writes recovery/conflict/recovery_jm.csv (per fit) and recovery_jm_sets.csv.
## From newsampler/recovery_jm.R, without its comparison with the legacy sampler's studies.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work")); OUT <- file.path(ROOT, "recovery", "conflict")
map <- read.csv(file.path(WORK, "rec_map.csv"), stringsAsFactors = FALSE)
rd <- function(f) { x <- readLines(f); k <- sub(":.*$", "", x); v <- sub("^[^:]*: ", "", x)
  keep <- !duplicated(k, fromLast = TRUE); setNames(v[keep], k[keep]) }
map$DIC <- NA_real_; map$rule <- NA_real_; map$met <- NA; map$meanD <- NA_real_
for (i in seq_len(nrow(map))) {
  lab <- if (file.exists(file.path(WORK, "summ", paste0(map$spec[i], "_jmhrec1b.rds")))) "_jmhrec1b" else "_jmhrec1"   # refit of record
  fs <- file.path(WORK, "summ", paste0(map$spec[i], lab, ".rds")); ft <- file.path(WORK, "fits", paste0(map$spec[i], lab, ".txt"))
  if (file.exists(fs)) { s <- readRDS(fs); if (is.data.frame(s$dic)) { map$DIC[i] <- s$dic$DIC[1]; map$meanD[i] <- s$dic$meanD[1] } }
  if (file.exists(ft)) { r <- rd(ft); g <- function(k) if (k %in% names(r)) r[[k]] else NA
    map$rule[i] <- suppressWarnings(as.numeric(g("stop_rule_value"))); map$met[i] <- g("stop_rule_met") %in% "TRUE" }
}
write_out(map[order(map$spec), ], file.path(OUT, "recovery_jm.csv"))
sets <- do.call(rbind, lapply(split(map, map$data_file), function(x) {
  done <- !anyNA(x$DIC); w <- if (done) which.min(x$DIC) else NA
  data.frame(study = x$study[1], gen = x$gen[1], generator = x$generator[1], replicate = x$replicate[1],
             n_cand = nrow(x), n_done = sum(!is.na(x$DIC)), n_met = sum(x$met %in% TRUE), true_met = isTRUE(x$met[x$true_model]),
             winner = if (done) x$candidate[w] else NA, recovered = if (done) x$true_model[w] else NA,
             gap = if (done) x$DIC[x$true_model] - x$DIC[w] else NA)
}))
sets <- sets[order(sets$study, sets$gen, sets$replicate), ]
write_out(sets, file.path(OUT, "recovery_jm_sets.csv"))
cat("fits with DIC:", sum(!is.na(map$DIC)), "of", nrow(map), "; stop rule met:", sum(map$met %in% TRUE), "\n")
for (s in unique(sets$study)) { y <- sets[sets$study == s & !is.na(sets$recovered), ]
  cat(sprintf("%-4s sets %2d: recovered %2d, largest gap %6.1f\n", s, nrow(y), sum(y$recovered), if (nrow(y)) max(y$gap) else NA)) }
