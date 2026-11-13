---
title: "Supplementary Materials: model fitting, model selection and graphics"
subtitle: "Racing to Control Decision Conflict"
author: "Andrew Heathcote, Marit Metz, Niek Stevenson & Dora Matzke"
---

These materials give the full detail of how the models in the paper were specified, fitted, checked, compared and plotted. They cover:

- the data and design (S1);
- the five model architectures' likelihoods (S2);
- the conflict mechanisms and their order restrictions (S3);
- the hierarchical model and priors (S4);
- the sampler (S5);
- convergence classes (S6);
- model selection (S7);
- posterior predictives and figures (S8);
- the model-recovery studies (S9);
- software and reproduction (S10).

Paths such as `code/07_figures/plot_heat_grid.R` refer to the reproducibility archive. Each section names the code that implements it.

Everything reported uses the new joint-move samplers described in S5 and the fits listed in `derived/tables/FITS.csv`. All predictives are at the subject level. All DICs are conditional DICs with the best-draw plug-in (S7).
