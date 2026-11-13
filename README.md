# ChoiceConflict: reproducibility archive

Heathcote, A., Metz, M., Stevenson, N., & Matzke, D. *Racing to Control Decision Conflict.*

This archive accompanies the paper's evidence-accumulation models (LBA, RDM, LNR, WDM, DDM) of simple choice and of the Simon, Stroop and attention-network (ANT) conflict tasks (data: Eisenberg et al., 2019). It holds:

- the supplement (`doc/`);
- the code that reproduces every fit, table and figure (`code/`);
- the EMC2/DMC2 builds that made the fits (`builds/`);
- the data, specs, fits, predictives and plot-ready tables (`inputs/`, `fits/`, `derived/`);
- the figures (`figures/`).

## Records

| record | holds |
|---|---|
| GitHub repository, https://github.com/andrewheathcote/ChoiceConflict, and its Zenodo software record (DOI: SOFTWARE_RECORD_DOI) | `code/`, `doc/`, `figures/`, `recovery/`, `builds/BUILDS.csv`, `MD5SUMS` and `install_builds.sh` |
| Zenodo data record (DOI: DATA_RECORD_DOI) | `inputs/`, the 41 fits behind the delta, CAF and defective-CDF figures (`derived/tables/figure_fits.csv`), their `derived/pp_raw/`, the mechanism knockouts (`derived/pp_null/`), `derived/tables/` and `derived/records/`, and the eleven EMC2/DMC2 source tarballs `builds/EMC2_<sha>.tar.gz` (`builds.tar`), as tar bundles |

The all-model figures (DIC, misfit and convergence grids) are drawn from the plot-ready tables, so their fits are not needed to reproduce them. The other 392 fits (~126 GB) are available from the authors on request.

To assemble a working archive: clone the repository and unpack the data record's bundles, including `builds.tar`, into the repository root (each bundle holds archive-relative paths).

## Layout

```
README.md  LICENSE (GPL-3)  LICENSE-CC-BY-4.0.txt  CITATION.cff  .zenodo.json
doc/       supplement.docx (the supplement), src/ (its Markdown + LaTeX source), make_tables.R, build.sh
code/      config.R (paths, builds), run_stage.R (one entry point per stage), lib/ (shared functions)
           00_data  01_specs  02_fit  03_classify  04_select  05_floors  06_predict  07_figures  08_recovery
builds/    EMC2_<sha>.tar.gz (11 builds), BUILDS.csv (which fits each made), install_builds.sh
inputs/    data/<task>/dat.RData (+ Stroop dat_acc.RData), specs/<spec>.rds, specs_unrestricted_q6/, noise_floors.csv
fits/      <task>/<arch>/<label>.RData   (local archive: all 433 = 353 conflict display fits, 24 simple-choice fits and 56 superseded runs used only by 03_classify/dic_tolerance.R; online: derived/tables/figure_fits.csv)
derived/   tables/ (plot-ready tables; FITS.csv lists every fit), pp_subj/ (subject-level predictive statistics),
           pp_raw/ (simulated trials for the EMC2-drawn figures), pp_null/ (the same with one mechanism removed: the
           paper's Figure 12), records/ (fit records), summ_refit/
figures/   the paper's and the supplement's figures (PDF); figures/README.md maps each file to its figure number
recovery/  conflict/ and choice/: design notes and results of the model-recovery studies (no fits, no simulated data)
```

## Install

1. R ≥ 4.4 with a C++17 compiler, plus the CRAN packages EMC2 depends on (Rcpp, RcppArmadillo, mvtnorm, coda, MASS, loo and others listed in each build's DESCRIPTION), and `parallel`.
2. Install the builds, each into its own library: unpack the data record's `builds.tar` into the archive root (it fills `builds/`), then `sh builds/install_builds.sh`, which checks them against `builds/MD5SUMS` and writes `lib/<sha8>/`.
   - `code/config.R` loads, for each stage, the build that made its objects, and checks the commit.
   - To install only what the figures need: `sh builds/install_builds.sh 3c1bcaf1 4f1fa123 f2f58985`.
3. Run from the archive root, or set `CC_ROOT` to it. Every path is relative to `CC_ROOT`; the other variables are listed at the top of `code/config.R`.

## Reproduce

| stage | command | time | needs |
|---|---|---|---|
| data | `Rscript code/run_stage.R data` | 20 min | the Eisenberg et al. trial files in `inputs/raw/` |
| classify | `Rscript code/run_stage.R classify` | seconds | `derived/tables/frozen/` |
| select | `Rscript code/run_stage.R select` | ~30 min | `fits/Choice/` |
| floors | `Rscript code/run_stage.R floors` | minutes | `recovery/`, `inputs/data/` |
| predict | `Rscript code/run_stage.R predict <cores>` | ~5 h on 8 cores (~11 GB RAM per process) | all displayed fits (local archive) |
| figures | `Rscript code/run_stage.R figures` | ~30 min | tables + `derived/pp_raw/` + `fits/Choice/` (choice CDFs) |
| one fit | `Rscript code/02_fit/fit.R <spec> <set> [cores] [chains]` | minutes (LNR) to days (ANT DDM) | `inputs/`, a build |
| recovery | `code/08_recovery/*/` (simulate, fit with sets jmh-rec1 / jmh-rec5 / jmdrc, summarise, score) | days | the generators' fits |

**Exact versus statistical reproduction.**

- *Exact:* with the online record alone, every table and figure reproduces. The all-model figures come from the plot-ready tables in `derived/tables/`; the delta/CAF figures and the choice CDFs come from `derived/pp_raw/`, `derived/pp_null/` and `fits/`.
- *Statistical:* a refit (`02_fit/fit.R`) starts from the original seed, but is not draw for draw: parallel chains, a different BLAS and the original runs' resumes after interruptions all change the random path. The fits ran on NCI Gadi, AWS and a laptop; the code here runs on any one machine.

**Fit sets.** `code/02_fit/fit_options.csv` lists every combination of build and sampler settings that made an archived fit, e.g. `jm`, `jmh-w2`, `jmd1`, except Niek Stevenson's 60 hard-case refits (set `niek`): those were made with his own driver on development builds, three of them on states not in the repository, so they are archived as fitted and are not refittable with `fit.R`. Each fit's set, build and chain count are in `derived/tables/FITS.csv`.

## Verify

What was checked when the archive was built (2026-11-07):

- **Fits.**
  - All 433 archived fits' md5s are in `derived/tables/MD5SUMS_fits.txt`.
  - The 240 fits that had a hash recorded at their source (Gadi, S3) match it.
- **Code against the original outputs.** Run from this archive alone, the ported code reproduced:
  - the classification table `display_jm.csv` (identical);
  - `choice_jm.csv` and the choice recovery scores (identical);
  - the delta/CAF floors (identical);
  - the convergence grids, the choice recovery figure and the study 5 recovery figure (0 differing pixels);
  - the DIC grids (only the cells changed by the Q6 refits differ);
  - all five analysed data files (Choice, Simon, Stroop, Stroop accuracy-coded, ANT) from the raw files (`code/00_data/make_data.R`; identical). This needed one fix to the cleaning script: its Stroop section reused the Simon section's excluded participants (s218, s375, s550). The analysed Stroop data, which every Stroop fit used, never excluded them, and the fixed script reproduces that data exactly.
- **DIC.** Recomputed from the archived fits, each on its own build (`compare(..., use_best_fit = TRUE)`), for the smallest fit of every conflict build and all 16 simple-choice fits: all 25 match the tables to within 0.001.
- **Figures.** Re-rendered from this archive (2026-11-08) and compared with the previous figures:
  - unchanged: choice, choice recovery, convergence grids, split6, models6 ANT and Stroop, best5 Stroop;
  - changed by the subject-level predictives, which replaced the earlier population-level ones: the delta and CAF grids, and three best5 picks (ANT CAF row WDM VcSc → BiVc; Simon WDM BiSc → Bi in the delta row and → BiVc in the CAF row);
  - changed by the restricted Simon LNR refits (S3): grid_dic_Simon and models6_Simon (the LNR MciSci panel); grid_dic_ANT and grid_dic_Stroop move slightly through the shared colour scale.
- **Online record alone.** The data record's bundles (`local_only/tools/stage_online.sh`), unpacked into an empty folder with only `code/` and `recovery/` added, render all 21 figures pixel-identical (0 differing pixels) to those in `figures/`.
- **Recovery tables.** The legacy-sampler comparison columns were removed.

## Licences

- Code: GPL-3 or later, as EMC2 (`LICENSE`).
- Data, tables, figures and text: CC-BY 4.0 (`LICENSE-CC-BY-4.0.txt`). The data derive from Eisenberg et al. (2019), released CC-BY 4.0; cite them.
