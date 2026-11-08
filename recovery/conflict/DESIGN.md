# Conflict-task model recovery: design, code and results

Neither the fits nor the simulated data sets are archived. The archive keeps:
- the rationale (this file);
- the code that simulates, fits and scores the studies (`code/08_recovery/conflict/`);
- the per-fit and per-set results (the CSV files here);
- the figures (`figures/recovery_jm_fig.pdf`, `figures/recovery_rec5_fig.pdf`).

## Why

Model selection by DIC needs a scale: how large must a DIC difference be before it means something? Each study simulates data from a known model fitted to the real data, then fits a set of candidates exactly as the real data were fitted, and asks two things: does the generating model win, and by how much does a wrong winner win? The largest wrong-winner margin over all studies is the **DIC floor** (149.4, rounded up to 150; `code/05_floors/recovery_floor.R`). The direction of the errors is also used: more complex wrong winners won by at most 40 (studies 1–4) and 28 (study 5), which gives the "=" rule's threshold of 50 for simpler models (supplement S7).

## How the data are simulated

The same procedure is used in every study:
- **Generator.** A model fitted to the real data. Studies 1–4 used earlier fits of the Simon models; study 5 uses the displayed fits.
- **Subject parameters.** For each real participant, draw a parameter vector from MVN(posterior mean of the population mean μ, posterior mean of the population covariance Σ).
- **Trials.** Simulate that participant's real trials, with the real design and trial counts, using `make_data(shrink2bound = TRUE)`. Out-of-bound parameters are clamped to the bound.
- **Replicates.** 5 per generator.

Each candidate is then fitted with the real model's specification (design, order restrictions, prior, chain count), with only the data swapped, and with the convergence campaign's sampler settings: fit sets `jmh-rec1` and `jmh-rec5` in `code/02_fit/fit_options.csv`, on build e8bf8344. The fits are scored by conditional DIC with the best-draw plug-in, exactly as the real fits.

## Studies 1–4 (Simon; 340 fits; `recovery_jm.csv`, `recovery_jm_sets.csv`)

| study | generators (Simon) | candidates | data sets |
|---|---|---|---|
| 1 LNR | MiSi, MiSci, MciSi, MciSci | the generators + McSc (no irrelevant-match terms: a control of the wrong shape) | 20 |
| 2 RDM | ViSi, ViSci, VciSi, VciSci | the generators + BcSci | 20 |
| 3 LBA | ViSi, ViSci, VciSi, VciSci | the generators + BcSci | 20 |
| 4 cross-architecture | RDM ViSi, DDM BiVc | RDM ViSi, RDM VciSci, DDM BiVc, DDM BciSc | 10 |

**Results:**
- Recovered: LNR 10/20, RDM 6/20, LBA 7/20, cross-architecture 10/10.
- Largest wrong-winner margins: 149.4, 81.6, 86.8 and 0.
- Every LNR and RDM miss is a simpler model winning. The LBA has three over-complex winners, with margins of at most 40.
- Cross-architecture: the generating architecture won every set:
  - by 955–1,288 (DDM data) and 2,105–2,424 (RDM data) over the other architecture's best;
  - by 29–124 over the nearest same-architecture model.

To rerun:
1. `studies1to4_specs.R <simdir>`
2. fit every spec with set `jmh-rec1`
3. `02_fit/summarise_fit.R` with `CC_DERIVED=<work>`
4. `score_studies1to4.R`

## Study 5: one mechanism or two? (280 fits; `recovery_rec5*.csv`)

For seven cells where the real data's best model has one mechanism or a close rival with two, the one-mechanism model (G1) and the two-mechanism model (G2) both generate data, and four candidates are fitted to each set: G1, G2 and two others.

| cell | task | architecture | G1 | G2 | other candidates |
|---|---|---|---|---|---|
| c1 | Simon | RDM | Bi | ViSi | BiSi, VcSi |
| c2 | Simon | LBA | Bi | Bci | BiSi, ViSi |
| c3 | Stroop | LNR | Mi | MiSi | Mci, MiSc |
| c4 | Stroop | LBA | Vi | Vci | BcVi, ViSi |
| c5 | ANT | RDM | Bi | BiVc | BiVi, ViSi |
| c6 | ANT | LBA | Bi | ViSi | BiSi, BiVi |
| c7 | Stroop | RDM | Vi | ViSi | Vci, BiSi |

**Results:**
- The one-mechanism generator was recovered in 32/35 sets, the two-mechanism generator in 21/35.
- Wrong winners' margins were at most 28 on G1 data and 53 on G2 data.
- For each cell, `recovery_rec5_cells.csv` gives false alarms, misses, and the real data's DIC difference against its simulated distributions.
- In Simon RDM (c1) and Simon LBA (c2) the real data's one-mechanism model is rejected, so those cells' one-mechanism models are never marked "=" (supplement S7).
- `dic_variants_rec5.csv` (`dic_variants.R`) shows why the best-draw DIC is used:
  - the standard DIC recovered 12/35 one-mechanism and 19/35 two-mechanism sets;
  - its worst wrong winner won by 8,961.
- RDM fits were made on e8bf8344, before the RDM likelihood fix. They were screened (`screen_old.rds` / `screen_fixed.rds` in the original project) and none was affected.

To rerun:
1. `rec5_simulate.R`
2. fit with set `jmh-rec5`
3. summarise
4. `score_rec5.R`, `dic_variants.R`, `recovery_figs.R`
