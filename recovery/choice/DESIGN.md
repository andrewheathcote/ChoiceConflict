# Simple-choice model recovery: design, code and results

## Why

On the simple-choice data the DDM and the LBA (mean rate free) are close: best-draw DIC favours the DDM by 93 and BPIC
favours the LBA by 51. Two questions follow:
- **Recoverability.** Can DIC and BPIC tell the two architectures apart at this data size, or could they mimic each other?
- **Size of a recoverable difference.** No DIC floor has been calibrated for choice; how large are the generating model's
  winning margins?

## Design

- **Generators.** The full DDM and LBA simple-choice fits (`Choice_DDM_v0_jmd1`, `Choice_LBA_jmd1x`, the LBA with its
  mean rate free), on their own build (4f1fa123).
- **Simulation.** Subject parameters drawn from MVN(posterior mean μ, posterior mean Σ), then `make_data(shrink2bound = TRUE)`
  on the 517 real participants and their trials. 5 replicates per generator, seed 20261104 + 1000·r + generator index
  (so the DDM data sets are those of the first study, below).
- **Candidates.** The best variant of each family: DDM_v0, LBA, LNR, RDM and WDM_v0 (race mean rate free). Each uses the
  real fit's spec, with only the data swapped. 50 fits.
- **Fitting.** Fit set `jmdrc2`: the simple-choice settings (at most 120 sample checks for the LBA and RDM, 40 for the
  others), on one 192-core AWS machine (run jmh-rc2, 2026-11-09; scripts, log and fit records in
  `local_only/historical/aws/jmh-rc2/`).
- **Scoring.** Best-draw conditional DIC and BPIC on 4f1fa123. `recch_dic.csv` was computed on the fitting machine as each
  fit finished (`rc2_dic.R`: the same `compare()` call as `dic.R`, plus the fit's sample draws and its true stop-rule value,
  `rule_value`); `score.R` scores it and draws `figures/choice_recovery.pdf`.

Code: `code/08_recovery/choice/`. Neither the fits nor the simulated data are archived.

## Results (`recch_dic.csv`, `recch_score.csv`, `figures/choice_recovery.pdf`)

- **Convergence.** 46 of 50 fits met the stop rule. Four race fits reached the 120-check cap without meeting it and are
  scored as fitted: LBA on DDM data r5 (rule 1.128), LBA on LBA data r2 (1.151) and r4 (1.122), RDM on LBA data r3 (1.175).
- **Recovery.** Both criteria recovered the generator in all 10 data sets; the runner-up was always the other of the DDM and LBA.
- **DDM − LBA difference** (LBA minus DDM; positive favours the DDM), equal here to the winning margin:

| data | DIC | BPIC |
|---|---|---|
| DDM-generated (5 sets) | 500 to 670 | 315 to 529 |
| LBA-generated (5 sets) | −679 to −885 | −965 to −1,307 |
| real data | 93 | −51 |

The real data's DDM–LBA differences are far smaller than any the study produced, and the two criteria disagree on their sign.

## Superseded first study (2026-11-04)

The first choice recovery study used the race models with the mean rate fixed (`Choice_LBA_v1_jmd1` as generator; LBA_v1
and RDM_v1 as candidates; fit set `jmdrc`, on Gadi) and also scored WAIC, PSIS-LOO and stacking. Every criterion recovered
the generator in all 10 sets. Its code and results are in `local_only/historical/choice_loo/` and `work/before_rc2/`.
WAIC and PSIS-LOO were dropped on 2026-11-10.
