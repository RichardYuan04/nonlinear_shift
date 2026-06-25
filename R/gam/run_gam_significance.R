# ---------------------------------------------------------------------------
# Region-wise GAM significance testing across the 360 HCP-MMP parcels.
#
# For each parcel we fit a generalized additive model
#     Activation ~ s(BIS) + Gender + Age
# and test whether the smooth term of task performance (BIS) explains
# additional deviance via an ANOVA-based likelihood ratio test against the
# nested linear model. p-values are FDR-corrected across all parcels.
#
# A second block fits a factor-smooth interaction to check whether the
# performance smooth interacts with the covariates (Age / Gender).
#
# Helper functions (gam.fit.smooth, gam.factorsmooth.interaction, ...) are
# defined in R/gam/gam_functions.R.
# ---------------------------------------------------------------------------

library(mgcv)
library(dplyr)
library(stringr)
library(gratia)

source("R/config.R")
source("R/gam/gam_functions.R")

# Model configuration ------------------------------------------------------
covariates <- c("Gender", "Age")
smooth_var <- "BIS"
int_var    <- "Age"   # covariate tested for a factor-smooth interaction
knots      <- 3

# ROI columns are the parcel activation columns; in the cleaned data they run
# from column 10 to the last column (360 parcels).
roi_from   <- 10

# Datasets to analyse. Place the corresponding CSVs in DATA_DIR (see
# data/README.md). Keys are used to name the output files.
dataset_files <- list(
  relational_relational = "relational_relational_cleaned.csv",
  WM_Body_2bk           = "WM_Body_2bk_cleaned.csv",
  WM_Face_2bk           = "WM_Face_2bk_cleaned.csv",
  WM_Place_2bk          = "WM_Place_2bk_cleaned.csv",
  WM_Tool_2bk           = "WM_Tool_2bk_cleaned.csv",
  math                  = "math_cleaned.csv",
  story                 = "story_cleaned.csv",
  ET_Shape              = "ET_Shape_cleaned.csv",
  ET_Face               = "ET_Face_cleaned.csv"
)

load_dataset <- function(file) read.csv(data_path(file))

# ===========================================================================
# Block 1: BIS x covariate (factor-smooth) interaction
# ===========================================================================
for (dataset_name in names(dataset_files)) {
  dataset     <- load_dataset(dataset_files[[dataset_name]])
  roi_columns <- colnames(dataset)[roi_from:ncol(dataset)]

  results        <- data.frame()
  gam.int.F      <- c()
  gam.int.pvalue <- c()

  for (roi in roi_columns) {
    result <- gam.factorsmooth.interaction(
      dataset, region = roi,
      smooth_var = smooth_var, int_var = int_var, covariates = covariates,
      knots = knots, set_fx = FALSE
    )
    results        <- rbind(results, result)
    gam.int.F      <- c(gam.int.F, result[1, "gam.int.F"])
    gam.int.pvalue <- c(gam.int.pvalue, result[1, "gam.int.pvalue"])
  }

  results$fdr_adjusted_pvalue <- p.adjust(gam.int.pvalue, method = "fdr")
  write.csv(results, file = result_path(paste0(dataset_name, "_BIS_age_results.csv")),
            row.names = FALSE)
}

# ===========================================================================
# Block 2: significance of the performance smooth term (main analysis)
# ===========================================================================
gam_datasets <- c("WM_Face_2bk", "WM_Place_2bk")

for (dataset_name in gam_datasets) {
  dataset     <- load_dataset(dataset_files[[dataset_name]])
  roi_columns <- colnames(dataset)[roi_from:ncol(dataset)]

  results        <- data.frame()
  var_p_values   <- c()
  anova_p_values <- c()

  for (roi in roi_columns) {
    result <- gam.fit.smooth(
      dataset, region = roi,
      smooth_var = smooth_var, covariates = covariates,
      knots = knots, set_fx = FALSE, stats_only = TRUE
    )
    results        <- rbind(results, result)
    var_p_values   <- c(var_p_values, result[1, "gam.smooth.pvalue"])
    anova_p_values <- c(anova_p_values, result[1, "anova.smooth.pvalue"])
  }

  results$var_fdr_adjusted_pvalue   <- p.adjust(var_p_values, method = "fdr")
  results$anova_fdr_adjusted_pvalue <- p.adjust(anova_p_values, method = "fdr")
  write.csv(results, file = result_path(paste0(dataset_name, "_gam_results.csv")),
            row.names = FALSE)
}
