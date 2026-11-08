# Simple-choice model recovery: design, code and results

## Why

On the simple-choice data the DDM is the best model by 147 DIC over the LBA. Two questions follow:
- **Recoverability.** Is the DDM–LBA difference recoverable, or could the two architectures mimic each other at this data size?
- **Agreement between criteria.** Do DIC, WAIC and PSIS-LOO agree? No DIC floor has been calibrated for choice.

## Design

- **Generators.** The full DDM and LBA simple-choice fits (`Choice_DDM_v0_jmd1`, `Choice_LBA_v1_jmd1`), on their own build (4f1fa123).
- **Simulation.** Subject parameters drawn from MVN(posterior mean μ, posterior mean Σ), then `make_data(shrink2bound = TRUE)` on the 517 real participants and their trials. 5 replicates per generator.
- **Candidates.** The best variant of each family: DDM_v0, LBA_v1, LNR, RDM_v1 and WDM_v0. Each uses the real fit's spec, with only the data swapped. 50 fits.
- **Fitting.** Fit set `jmdrc`: the simple-choice settings, on Gadi.
- **Scoring.** Best-draw conditional DIC on 4f1fa123 (`dic.R`); WAIC, PSIS-LOO, pseudo-BMA+ and stacking on DMC2 main f2f58985 (`loo.R`, with `use_best_fit = TRUE` so its DIC matches); `score.R`, `score_recch_ic.R`.

Code: `code/08_recovery/choice/`. Neither the fits nor the simulated data are archived.

## Results (`recch_dic.csv`, `recch_score.csv`, `recch_ic_score.csv`, `loo/`, `figures/choice_recovery_ic.pdf`)

- **Recovery.** Every criterion recovered the generator in all 10 data sets.
- **DDM − LBA difference** (positive favours the DDM):

| data | DIC | WAIC | PSIS-LOO |
|---|---|---|---|
| DDM-generated (5 sets) | 648 to 834 | 872 to 1,095 | 861 to 1,076 |
| LBA-generated (5 sets) | −936 to −1,294 | −799 to −1,198 | −864 to −1,236 |
| real data | 147 | 580 | 502 |

The real data's DDM advantage is smaller than in either generator's simulations, but positive under every criterion.
