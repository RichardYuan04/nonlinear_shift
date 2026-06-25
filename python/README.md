# Python notebooks

Downstream analyses that consume the per-parcel derivative / GAM results
produced by the R pipeline (see `../R`). Install dependencies with
`pip install -r ../requirements.txt`.

The notebooks expect their input CSVs (per-parcel derivative results, GMM
outputs, the S-A ranking table, etc.) in the repository `data/` directory;
adjust the `data/...` paths at the top of each notebook if you keep them
elsewhere. Subject-level inputs must be obtained from HCP (see
`../data/README.md`). All cell outputs have been cleared.

| Notebook | Purpose |
|----------|---------|
| `01_gmm_linearity_filter.ipynb` | Fit a Gaussian Mixture Model to the second-derivative distribution and drop near-linear (≈0) parcels, keeping only nonlinear regions. |
| `02_sa_rank_correlation.ipynb` | Spearman correlation between parcel nonlinearity (mean 2nd derivative) and S-A rank, by condition. |
| `03_permutation_tests.ipynb` | Label-shuffling permutation tests comparing S-A ranks and concave/convex counts across conditions. |
| `04_inflection_points.ipynb` | Locate inflection (turning) points of concave vs. convex curves and compare them. |
| `05_ttest_anova.ipynb` | Behavioral and inflection-point group comparisons (t-tests, ANOVA, normality / variance checks). |
| `06_ml_decoding.ipynb` | Decode performance (BIS) from parcel activation via LASSO, logistic regression, and SVR. |
| `07_random_split_replication.ipynb` | Replicate the S-A correlation on two random split-half subsamples. |
