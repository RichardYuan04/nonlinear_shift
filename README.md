# Nonlinear shift along the sensorimotor-association axis in brain responses to task performance

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.neuroimage.2026.121978-blue)](https://doi.org/10.1016/j.neuroimage.2026.121978)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey)](http://creativecommons.org/licenses/by/4.0/)

> 📄 **Official code repository for our paper published in *NeuroImage* (2026, open access).**
> Yuan, Y., Zhang, B., Perkins, K., & Cao, F. (2026). *Nonlinear shift along the
> sensorimotor-association-axis in brain responses to task performance.*
> **NeuroImage, 334, 121978.** https://doi.org/10.1016/j.neuroimage.2026.121978
> See [Citation](#-citation) below.

This project investigates how cortical activation across 360 brain regions responds to working memory (WM) task performance under four visual stimulus categories (**faces, places, tools, and body parts**), using high-quality fMRI data from the Human Connectome Project (HCP). We characterize how the brain–performance relationship changes **non-linearly** and how this change is organized along the **sensorimotor-to-association (S-A) axis**.

## 🧠 Project Motivation

Cognitive neuroscience often assumes a *linear* relationship between brain activation and task performance, yet findings conflict across studies. We hypothesized that, across a full range of performance, the relationship may be **non-linear**, and that the form of this relationship varies **systematically along the S-A axis** — with higher-rank association regions showing greater concavity (an inverted-U) than lower-rank sensorimotor regions. We further asked whether this pattern is modulated by visual stimulus category.

## 🔑 Key Findings

- A **gradual shift along the S-A axis**: higher-rank regions show greater concavity (inverted-U) than lower-rank regions — but **only in the face and body-part conditions**.
- In the **place** condition, many regions instead show a **convex** (U-shaped) pattern; in **tool/place**, few high-order regions engage, so the S-A correlation is absent.
- **Inflection points** sit *above* average performance in concave regions and *below* average in convex regions.
- Reveals a **novel principled mapping from brain topology (S-A rank) to brain function** in response to task proficiency.

## 📊 Methodology

- **Data Source:** HCP task-fMRI data (MSM-All registered, HCP-MMP 1.0 parcellation, 360 parcels)
- **Behavioral Measure:** BIS (Balanced Integration Score) = standardized accuracy − standardized RT, robust to speed–accuracy trade-offs
- **Modeling:** Region-wise Generalized Additive Models — `Activation ~ s(BIS) + Gender + Age`
- **Significance:** ANOVA-based likelihood ratio test (`anova.gam`) with False Discovery Rate (FDR) correction (P_FDR < 0.05)
- **Pattern classification:** First- and second-order derivatives of the fitted curves, with a Gaussian Mixture Model (GMM) separating linear vs. non-linear regions → six response types (linear/concave/convex × increase/decrease)
- **S-A axis:** Sensorimotor-Association ranking following Sydnor et al. (2021)

## 📁 Project Structure

```
├── BrainMap/
│   └── region_type.R        # Brain region visualization using ggseg and ggplot2
│
├── GAM/
│   ├── aic.R                # AIC-based model selection
│   ├── fdr_correction.R     # ANOVA-based testing and FDR correction for GAM models
│   ├── gam_functions.R      # Helper functions for GAM fitting, derivative calculation, etc.
│   └── plot_gam.R           # Visualization of fitted GAM curves and scatter plots
```

## 📌 Dependencies

- R (≥ 4.0)
- `mgcv`
- `ggplot2`
- `ggseg`
- `ggpubr`
- `dplyr`, `purrr`, `tidyr`
- Optional: `fdrtool`, `boot` for significance testing

## 📈 Example Output

- Cortical maps showing region-wise effect size and derivative-based region-type classifications
- GAM fit results
- GAM-fitted activation curves by BIS

## 📂 Data Availability

Data are from the **Human Connectome Project (HCP) Young Adult dataset (S1200 release)**, publicly available under the HCP Open Access and Restricted Access Data Use Terms. This repository contains the analysis code; it does not redistribute HCP data.

## 📄 Citation

If you use any part of the code, methods, or findings, please cite:

> Yuan, Y., Zhang, B., Perkins, K., & Cao, F. (2026). Nonlinear shift along the
> sensorimotor-association-axis in brain responses to task performance.
> *NeuroImage, 334*, 121978. https://doi.org/10.1016/j.neuroimage.2026.121978

**BibTeX:**

```bibtex
@article{yuan2026nonlinear,
  title     = {Nonlinear shift along the sensorimotor-association-axis in brain responses to task performance},
  author    = {Yuan, Yuqi and Zhang, Bohan and Perkins, Kyle and Cao, Fan},
  journal   = {NeuroImage},
  volume    = {334},
  pages     = {121978},
  year      = {2026},
  publisher = {Elsevier},
  doi       = {10.1016/j.neuroimage.2026.121978}
}
```

## 🙏 Acknowledgments

- Human Connectome Project (HCP) for data access
- HCP-MMP 1.0 parcellation atlas (Glasser et al., 2016)
- R packages: `mgcv`, `ggseg`, `ggplot2`, etc.

## 📬 Contact

- **Yuqi Yuan** (first author) — u3597424@connect.hku.hk — [@RichardYuan04](https://github.com/RichardYuan04)
- **Fan Cao** (corresponding author) — fancao@hku.hk
- **Bohan Zhang** — zhangbohan@chd.edu.cn

---

This repository accompanies our *NeuroImage* (2026) paper. Released under CC BY 4.0 — please cite as above if you use the code, methods, or results.
