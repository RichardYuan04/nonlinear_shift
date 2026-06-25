# Nonlinear shift along the sensorimotor–association axis

Analysis code for:

> Yuan, Y., Zhang, B., Perkins, K., & Cao, F. (2026). *Nonlinear shift along
> the sensorimotor–association axis in brain responses to task performance.*
> **NeuroImage**, 334, 121978. https://doi.org/10.1016/j.neuroimage.2026.121978

## Overview

Cognitive neuroscience often assumes a **linear** relationship between brain
activation and task performance, yet reported findings conflict. Using the
Human Connectome Project (HCP) working-memory task, we model activation in each
of the 360 HCP-MMP cortical parcels as a **nonlinear** function of individual
task performance and ask whether the *shape* of that function varies
systematically along the **sensorimotor–association (S-A) axis** and across
four visual stimulus categories (faces, places, tools, body parts).

Key findings: higher-rank (association) regions show greater **concavity**
(inverted-U) with performance — but only in the face and body-part conditions;
in the place condition many lower-rank regions show a **convex** (U-shaped)
pattern; and inflection points sit above average performance for concave
regions and below it for convex regions.

## Analysis pipeline

| Stage | What it does | Code |
|-------|--------------|------|
| 1. GAM significance | Region-wise `Activation ~ s(BIS) + Gender + Age`; ANOVA-based smooth test + FDR correction | `R/gam/run_gam_significance.R`, `R/gam/gam_functions.R` |
| 2. Model selection | AIC grid search over the smooth basis `k` and penalty order `m` | `R/gam/model_selection_aic.R` |
| 3. Derivatives | Bootstrap significance of the mean 1st/2nd derivative → direction & curvature per parcel | `R/derivatives/bootstrap_derivatives.R` |
| 3b. Robustness | Nonparametric derivative cross-check (categorical regression splines / local polynomial) | `R/derivatives/crs_spline_derivatives.R`, `R/derivatives/local_poly_derivatives.R`, `R/derivatives/boot_der.R` |
| 4. Curve / brain maps | Fitted GAM curves; cortical surface maps of response type | `R/gam/plot_gam.R`, `R/brainmap/region_type_brainmap.R` |
| 5. Linearity filter | Gaussian Mixture Model to drop near-linear (≈0 second-derivative) parcels | `python/01_gmm_linearity_filter.ipynb` |
| 6. S-A correlation | Spearman correlation between nonlinearity and S-A rank | `python/02_sa_rank_correlation.ipynb` |
| 7. Permutation tests | Label-shuffling tests comparing S-A ranks / curvature counts across conditions | `python/03_permutation_tests.ipynb` |
| 8. Inflection points | Locate and compare turning points of concave vs. convex curves | `python/04_inflection_points.ipynb` |
| 9. t-test / ANOVA | Behavioral and inflection-point group comparisons | `python/05_ttest_anova.ipynb` |
| 10. ML decoding | LASSO / logistic / SVR decoding of performance from activation | `python/06_ml_decoding.ipynb` |
| 11. Replication | Re-run the S-A correlation on two random split-half subsamples | `python/07_random_split_replication.ipynb` |

## Repository layout

```
.
├── R/                         # GAM modeling, derivatives, brain maps (run from repo root)
│   ├── config.R               # DATA_DIR / RESULTS_DIR convention
│   ├── gam/                   # GAM fitting, significance, model selection, plotting
│   ├── derivatives/           # bootstrap & nonparametric derivative tests
│   └── brainmap/              # ggseg cortical surface maps
├── python/                    # Jupyter notebooks for the downstream analyses
├── data/                      # data template + how to obtain the HCP data (inputs git-ignored)
└── results/                   # generated outputs (git-ignored)
```

## Data availability

Data are from the **HCP Young Adult dataset (S1200)** and are **not** included
here, per the HCP Data Use Terms. See [`data/README.md`](data/README.md) for
how to obtain the data and the expected input schema; `data/template_WM_2bk.csv`
shows the column layout with synthetic rows.

## Requirements

- **R** (≥ 4.0): `mgcv`, `gratia`, `dplyr`, `tidyr`, `stringr`, `boot`,
  `ggplot2`, `ggseg`, `ggsegGlasser`, `crs`, `parabar`, `doParabar`, `foreach`
- **Python** (≥ 3.9): see [`requirements.txt`](requirements.txt)

R scripts are designed to be run from the repository root (e.g.
`Rscript R/gam/run_gam_significance.R`). Inputs are read from `data/` and
outputs written to `results/`.

## Citation

If you use this code, please cite the paper above.

## Contact

- **Yuqi Yuan** — u3597424@connect.hku.hk — [RichardYuan04](https://github.com/RichardYuan04)
- **Bohan Zhang** — zhangbohan@chd.edu.cn

## License

Code released under the [MIT License](LICENSE).
