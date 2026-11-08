## Conditional best-draw DIC for the 16 simple-choice fits (label jmd1), on their own build (4f1fa123). Read-only.
##   Rscript code/04_select/choice_dic.R      -> derived/tables/choice_dic.csv (then choice_tables.R -> choice_jm.csv)
## compare(..., BayesFactor = FALSE, WAIC = FALSE, type = "conditional", use_best_fit = TRUE) per fit, as for every fit in
## the project; Dmean is the deviance at the posterior mean (the figure's "D"). From newsampler/gadi/choice/choice_dic.R.
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
SHA <- use_build(BUILD[["choice"]])
FT <- read_tab("FITS.csv"); FT <- FT[FT$task == "Choice", ]
rows <- lapply(seq_len(nrow(FT)), function(i) {
  e <- load_fit(FT$label[i])
  ic <- compare(list(m = e), BayesFactor = FALSE, WAIC = FALSE, type = "conditional", use_best_fit = TRUE, print_summary = FALSE)
  cat(sprintf("%-16s DIC %.1f Dmean %.1f\n", FT$spec[i], ic$DIC, ic$Dmean))
  data.frame(spec = FT$spec[i], variant = sub("^Choice_", "", sub("_v[01]$", "", FT$spec[i])), label = FT$label[i],
             n_par = length(e[[1]]$par_names), ic, row.names = NULL, check.names = FALSE)
})
D <- do.call(rbind, rows); D$emc2_sha <- SHA
D <- D[order(D$DIC), ]; D$dDIC <- D$DIC - min(D$DIC)
write_out(D, tab_path("choice_dic.csv"))
print(D[, c("variant", "n_par", "DIC", "dDIC", "Dmean", "EffectiveN")], row.names = FALSE)
