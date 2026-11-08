## Q6 (2026-11-07): the eight Simon LNR models with a congruency term on the location m (Mc, McSc, McSci, McSi, Mci,
## MciSc, MciSci, MciSi) get the order restriction m_CIincong = exp(.) (incongruent trials finish later), as in ANT and
## Stroop; m_iMd stays restricted where present. The original fits estimated m_CIincong freely (S3).
##   Rscript code/01_specs/restrict_q6.R     inputs/specs_unrestricted_q6/ -> inputs/specs/ (the restricted specs)
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
for (v in c("Mc", "McSc", "McSci", "McSi", "Mci", "MciSc", "MciSci", "MciSi")) {
  s <- readRDS(file.path(INPUTS, "specs_unrestricted_q6", paste0("Simon_LNR_", v, ".rds")))
  stopifnot("m_CIincong" %in% s$stored_par_names)
  f <- if (is.null(s$pre)) character(0) else s$pre$func
  f["m_CIincong"] <- "exp"; s$pre <- list(func = f)
  s$restriction_note <- "m_CIincong=exp added 2026-11-07 (Q6); the stored legacy fit had it free"
  saveRDS(s, file.path(INPUTS, "specs", paste0("Simon_LNR_", v, ".rds")))
  cat(v, ":", paste(names(f), f, sep = "=", collapse = " "), "\n")
}
