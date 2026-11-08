## Conditional best-draw DIC for the 50 choice recovery fits (label jmdrc), on their own build (4f1fa123), exactly as
## 04_select/choice_dic.R scores the real fits. Read-only on the fits.
##   Rscript code/08_recovery/choice/dic.R      -> recovery/choice/recch_dic.csv
## From newsampler/gadi/choice_rec/recch_dic.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
SHA <- use_build(BUILD[["choice"]])
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work"))
map <- read.csv(file.path(WORK, "recch_map.csv"), stringsAsFactors = FALSE)
rows <- lapply(map$spec, function(s) {
  f <- file.path(WORK, "fits", paste0(s, "_jmdrc.RData")); txt <- readLines(file.path(WORK, "fits", paste0(s, "_jmdrc.txt")))
  ev <- new.env(); nm <- load(f, envir = ev); e <- ev[[nm[1]]]; rm(ev)
  ic <- compare(list(m = e), BayesFactor = FALSE, WAIC = FALSE, type = "conditional", use_best_fit = TRUE, print_summary = FALSE)
  data.frame(spec = s, n_par = length(e[[1]]$par_names), ic,
             stop_rule_met = sub("^stop_rule_met: ", "", grep("^stop_rule_met:", txt, value = TRUE)[1]),
             stop_rule_value = as.numeric(sub("^stop_rule_value: ", "", grep("^stop_rule_value:", txt, value = TRUE)[1])),
             row.names = NULL, check.names = FALSE)
})
D <- merge(map, do.call(rbind, rows), by = "spec"); D$emc2_sha <- SHA
write_out(D, file.path(ROOT, "recovery", "choice", "recch_dic.csv")); cat("done:", nrow(D), "of", nrow(map), "\n")
