## One entry point per stage. Run from anywhere:  Rscript code/run_stage.R <stage> [cores]
## Stages, in order (each needs the outputs of the ones before it; the archive ships every stage's outputs, so any stage
## can be rerun on its own):
##   data       00_data/make_data.R               raw Eisenberg et al. csv files -> inputs/data/<task>/dat.RData
##   classify   03_classify/{display_index, apply_q6, mech_key}.R   frozen classification inputs -> display_jm.csv, mech_key.csv
##   select     04_select/{choice_dic, choice_tables, threshold_compare}.R
##   floors     05_floors/{recovery_floor, noise_floor_agg}.R
##   predict    06_predict/pp_subj.R run + raw, build_adequacy.R (mean and median), pp_null.R (mechanism knockouts)
##   figures    07_figures/*.R + 08_recovery/*/figures                 -> figures/
## Not run by this script (hours to days of computing, per model): 01_specs (specs from stored fits), 02_fit/fit.R (one fit
## per call), 02_fit/summarise_fit.R, 03_classify/dic_tolerance.R (reads every conflict fit, so local archive only; ~10 min on
## 4 cores), and the recovery studies' simulate / fit / score scripts (08_recovery/*/).
CODE <- local({ d <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))); while (!file.exists(file.path(d, "config.R"))) d <- dirname(d); d })
a <- commandArgs(TRUE); stage <- a[1]; cores <- if (length(a) > 1) a[2] else "4"
run <- function(script, ..., env = character()) {
  cmd <- paste(c(env, "Rscript", shQuote(file.path(CODE, script)), ...), collapse = " ")
  cat("\n==>", cmd, "\n"); st <- system(cmd); if (st != 0) stop(script, " failed (exit ", st, ")")
}
switch(stage,
  data     = run("00_data/make_data.R"),
  classify = { run("03_classify/display_index.R"); run("03_classify/apply_q6.R"); run("03_classify/mech_key.R") },
  select   = { run("04_select/choice_dic.R"); run("04_select/choice_tables.R"); run("04_select/threshold_compare.R") },
  floors   = { run("05_floors/recovery_floor.R"); run("05_floors/noise_floor_agg.R") },
  predict  = { run("06_predict/pp_subj.R", "run", cores)
               run("06_predict/build_adequacy.R"); run("06_predict/build_adequacy.R", env = "PP_AGG=median")
               run("07_figures/raw_needed.R", cores)                      # pp_raw for every figure fit (figure_fits.csv)
               run("06_predict/pp_null.R", cores, rbind(rep(c("Simon_RDM_ViSi_jm", "Stroop_RDM_ViSi_jm", "ANT_RDM_ViSi_jmx"),
                   each = 2), c("v_iMd", "s_iMd"))) },                     # the paper's Figure 12 knockouts
  ## the paper's and the supplement's figures only (figures/README.md maps each file to its figure number)
  figures  = { run("07_figures/plot_heat_grid.R", "dic")                                         # paper 7-9, S1
               for (m in c("delta", "caf")) run("07_figures/plot_heat_grid.R", m, env = "PP_AGG=median")   # S2, S3
               run("07_figures/plot_lattice.R"); run("07_figures/plot_models6.R")                 # S4, S5
               run("07_figures/plot_best5_grid.R", "t")                                           # paper 10
               run("07_figures/plot_models6_mix.R", "Stroop", "ANT", "RDM")                       # paper 11
               run("07_figures/plot_mech_null_grid.R")                                            # paper 12
               run("07_figures/plot_split6_grid.R")                                               # paper 13
               for (m in c("RDM_BiSi", "RDM_VciSci", "RDM_BciSci", "LNR_MciSci")) run("07_figures/plot_split6.R", 100, m)   # S6
               run("07_figures/plot_effects.R")                                                   # paper 14
               run("07_figures/plot_choice.R", "plot"); run("07_figures/plot_choice_types.R")     # paper 4, S7; S8
               run("08_recovery/conflict/recovery_figs.R"); run("08_recovery/choice/score.R") },  # S9, S10; S11
  stop("unknown stage ", stage, ": data, classify, select, floors, predict, figures"))
