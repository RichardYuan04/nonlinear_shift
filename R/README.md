# R analysis scripts

GAM modeling, derivative significance testing, and cortical surface
visualization. Run every script **from the repository root** so the relative
`data/` and `results/` paths resolve, e.g.:

```bash
Rscript R/gam/run_gam_significance.R
```

`config.R` defines `DATA_DIR` (default `data/`) and `RESULTS_DIR` (default
`results/`); override them with the `NLS_DATA_DIR` / `NLS_RESULTS_DIR`
environment variables.

| File | Purpose |
|------|---------|
| `config.R` | Shared data/results path configuration and helpers. |
| `gam/gam_functions.R` | Core GAM helpers (mgcv): fit the performance smooth, partial R², ANOVA test, posterior smooths, derivatives. `source()` before use. |
| `gam/run_gam_significance.R` | Region-wise GAM significance + FDR across the 360 parcels (main analysis); plus a BIS × covariate interaction check. |
| `gam/model_selection_aic.R` | AIC grid search over smooth basis `k` and penalty order `m`. |
| `gam/plot_gam.R` | Plot a fitted activation–performance curve with 95% CI for one parcel. |
| `derivatives/bootstrap_derivatives.R` | Bootstrap test of the mean 1st/2nd derivative per parcel (direction & curvature). |
| `derivatives/boot_der.R` | Residual-bootstrap helper (`compute_mean_deriv`) following Racine (1997) / the `crs` package. |
| `derivatives/crs_spline_derivatives.R` | Robustness check: derivatives via categorical regression splines (`crs::npglpreg`). |
| `derivatives/local_poly_derivatives.R` | Exploratory local-polynomial derivative estimation and diagnostic plots. |
| `brainmap/region_type_brainmap.R` | ggseg cortical maps coloring each parcel by its response type. |
