## Conflict recovery studies 1-4: the candidates' specs for each simulated data set.
##   Rscript code/08_recovery/conflict/studies1to4_specs.R <simdir>
##     simdir  folder with the simulated Simon data sets sim_<gen>_r<rep>.RData of the four studies, one subfolder each:
##             lnr/, rdm/, lba/, xa/ (not archived; simulated as in rec5_simulate.R, from fits of the generating models:
##             subject parameters from MVN(posterior mean mu, posterior mean Sigma), make_data(shrink2bound = TRUE) on the
##             real Simon subjects and trials, 5 replicates per generator; see recovery/conflict/DESIGN.md)
## Studies (generators; candidates = the generators plus one control; 5 replicates each):
##   LNR  MiSi MiSci MciSi MciSci; + McSc        RDM  ViSi ViSci VciSi VciSci; + BcSci
##   LBA  ViSi ViSci VciSi VciSci; + BcSci       XA   RDM_ViSi, DDM_BiVc; candidates RDM_ViSi RDM_VciSci DDM_BiVc DDM_BciSc
## 100 + 100 + 100 + 40 = 340 fits. Each candidate is the Simon model's own spec (design, prior, chain count) with only the
## data swapped. Writes <work>/inputs/{specs,data/Simon} and <work>/rec_map.csv; fit each spec with set jmh-rec1
## (CC_INPUTS=<work>/inputs CC_REFITS=<work>/fits Rscript code/02_fit/fit.R <spec> jmh-rec1), then score_studies1to4.R.
## From newsampler/specs_rec/make_specs_rec.R. Spec names r<study>_g<k>_r<rep>_<candidate>: k indexes the generator
## (never a mechanism string, which can collide in letter case on case-insensitive file systems).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
source(file.path(CODE, "config.R"))
simdir <- commandArgs(TRUE)[1]; stopifnot(!is.na(simdir), dir.exists(simdir))
WORK <- Sys.getenv("CC_RECOVERY_WORK", file.path(ROOT, "recovery_work"))
for (d in c("inputs/specs", "inputs/data/Simon")) dir.create(file.path(WORK, d), recursive = TRUE, showWarnings = FALSE)
studies <- list(
  LNR = list(dir = "lnr", gen = c(g1 = "gMS", g2 = "gMSs", g3 = "gMmS", g4 = "gMmSs"),
             truth = c("LNR_MiSi", "LNR_MiSci", "LNR_MciSi", "LNR_MciSci"), cand = c("LNR_MiSi", "LNR_MiSci", "LNR_MciSi", "LNR_MciSci", "LNR_McSc")),
  RDM = list(dir = "rdm", gen = c(g1 = "g1_vs", g2 = "g2_vSs", g3 = "g3_Vvs", g4 = "g4_VvSs"),
             truth = c("RDM_ViSi", "RDM_ViSci", "RDM_VciSi", "RDM_VciSci"), cand = c("RDM_ViSi", "RDM_ViSci", "RDM_VciSi", "RDM_VciSci", "RDM_BcSci")),
  LBA = list(dir = "lba", gen = c(g1 = "gvs", g2 = "gvSs", g3 = "gVvs", g4 = "gVvSs"),
             truth = c("LBA_ViSi", "LBA_ViSci", "LBA_VciSi", "LBA_VciSci"), cand = c("LBA_ViSi", "LBA_ViSci", "LBA_VciSi", "LBA_VciSci", "LBA_BcSci")),
  XA  = list(dir = "xa", gen = c(g1 = "gRDM", g2 = "gDDM"),
             truth = c("RDM_ViSi", "DDM_BiVc"), cand = c("RDM_ViSi", "RDM_VciSci", "DDM_BiVc", "DDM_BciSc")))
map <- NULL
for (s in names(studies)) { st <- studies[[s]]
  for (r in 1:5) for (k in seq_along(st$gen)) {
    src <- file.path(simdir, st$dir, sprintf("sim_%s_r%d.RData", st$gen[[k]], r)); stopifnot(file.exists(src))
    df <- sprintf("r%s_%s_r%d.RData", s, names(st$gen)[k], r)
    file.copy(src, file.path(WORK, "inputs/data/Simon", df), overwrite = FALSE)
    for (cand in st$cand) {
      sp <- read_spec(paste0("Simon_", cand)); nm <- sprintf("r%s_%s_r%d_%s", s, names(st$gen)[k], r, cand)
      sp$data_file <- df
      sp$recovery <- list(study = s, generator = st$truth[k], replicate = r, candidate = cand, true_model = cand == st$truth[k])
      saveRDS(sp, file.path(WORK, "inputs/specs", paste0(nm, ".rds")))
      map <- rbind(map, data.frame(spec = nm, study = s, gen = names(st$gen)[k], generator = st$truth[k], replicate = r,
                                   candidate = cand, true_model = cand == st$truth[k], arch = sp$arch, chains = sp$n_chains, data_file = df))
    }
  }
}
stopifnot(!anyDuplicated(tolower(map$spec)))
write.csv(map, file.path(WORK, "rec_map.csv"), row.names = FALSE)
cat(nrow(map), "specs;", sum(map$arch != "DDM"), "race,", sum(map$arch == "DDM"), "DDM\n")
