# Data

This study is a secondary analysis of the **Human Connectome Project (HCP)
Young Adult dataset (S1200 release)**. Subject-level data are governed by the
[HCP Open Access and Restricted Access Data Use Terms](https://www.humanconnectome.org/study/hcp-young-adult/data-use-terms)
and **are not redistributed in this repository**. Only code and a synthetic
data template are provided here.

## How to obtain the data

1. Register for ConnectomeDB and accept the HCP Data Use Terms:
   https://db.humanconnectome.org
2. Download the WM task contrasts (2-back vs. fixation) and the unrestricted /
   restricted behavioral data for the S1200 release.
3. Reproduce the per-subject, per-parcel activation tables (HCP-MMP 1.0
   parcellation, 360 cortical parcels) and the behavioral measures, then place
   the resulting CSVs in this `data/` directory using the file names referenced
   by the scripts (e.g. `WM_Face_2bk_cleaned.csv`, `WM_Place_2bk_cleaned.csv`,
   `WM_Body_2bk_cleaned.csv`, `WM_Tool_2bk_cleaned.csv`).

All scripts read inputs from this directory and write outputs to `results/`
(both are git-ignored). The R scripts honor the `NLS_DATA_DIR` and
`NLS_RESULTS_DIR` environment variables if you prefer a different location.

## Input schema

`template_WM_2bk.csv` shows the expected column layout with **three synthetic
rows** (placeholder subject IDs, random values — not real data). Each cleaned
input CSV has 369 columns:

| Column        | Description                                                            |
|---------------|------------------------------------------------------------------------|
| `Subject`     | Subject identifier                                                     |
| `Gender`      | `M` / `F` (linear covariate)                                           |
| `Age`         | Age in years (linear covariate)                                        |
| `pooled_Acc`  | Accuracy pooled across conditions                                      |
| `<Cond>_Acc_z`| Standardized (z-scored) accuracy for the condition                     |
| `<Cond>_RT_z` | Standardized (z-scored) median RT for the condition                    |
| `BIS`         | Balanced Integration Score = standardized accuracy − standardized RT   |
| `<Cond>_ACC`  | Raw accuracy                                                           |
| `<Cond>_RT`   | Raw median reaction time (ms) on correct trials                        |
| `*_ROI` (×360)| Mean parameter estimate (activation) per HCP-MMP 1.0 cortical parcel   |

`BIS` is the smooth (predictor) variable; the 360 `*_ROI` columns are the
response variables, modeled one at a time. Parcel columns run from column 10 to
column 369.

## S-A axis ranking

The sensorimotor–association (S-A) axis ranking used in the correlation
analyses is derived following Sydnor et al. (2021), integrating 10
neurobiological maps. The ranking table (`region` → `final.rank`, 1 = most
sensorimotor, 180 = most association) is keyed by parcel name and is not
subject-level data.
