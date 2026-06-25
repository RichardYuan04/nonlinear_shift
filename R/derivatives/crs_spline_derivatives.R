# ---------------------------------------------------------------------------
# Robustness analysis: mean first/second derivatives via categorical
# regression splines (crs::npglpreg) instead of mgcv GAMs.
#
# For each parcel we select bandwidths by least-squares cross-validation, then
# use the residual-bootstrap procedure in boot_der.R (compute_mean_deriv) to
# obtain a two-sided p-value for the mean first and second derivative of the
# performance smooth. This provides a nonparametric cross-check of the GAM
# derivative results. Run from the repository root.
# ---------------------------------------------------------------------------

# Load required R packages
library("crs")
library("parabar")
library("doParabar")
library("foreach")

source("R/config.R")

# Load the compute_mean_deriv() helper
source("R/derivatives/boot_der.R")

# Read data (place the cleaned CSV in DATA_DIR; see data/README.md)
df <- read.csv(data_path("WM_Face_2bk_cleaned.csv"))
vars <- colnames(df)
Y_names <- vars[10:369]  # response (parcel) columns: columns 10..369 (adjust per dataset)

# Extract predictors
X <- df[c("BIS", "Age", "Gender")]
X$Age <- ordered(X$Age)      # treat Age as an ordered factor
X$Gender <- factor(X$Gender) # treat Gender as a factor

# Initialize the results data frame
results <- data.frame(
  label = character(),
  First_Mean = numeric(),
  First_P = numeric(),
  First_Sig = logical(),
  Second_Mean = numeric(),
  Second_P = numeric(),
  Second_Sig = logical(),
  stringsAsFactors = FALSE
)

# Loop over each response (parcel) variable
for (Y in Y_names) {
  cat(sprintf("Processing variable: %s\n", Y))

  # Build the model formula
  f <- formula(paste(Y, "~ BIS + ordered(Age) + factor(Gender)"))

  # Automatic bandwidth selection
  model_bw <- npglpreg(formula = f, data = df, bwtype = "auto", degree = 4, nmulti = 5, cv = "bandwidth")
  bws <- model_bw$bws

  # First derivative
  result_first <- compute_mean_deriv(
    Y = df[, Y],
    X = X,
    bws = bws,
    degree = 4,
    gradient.vec = 1,  # first derivative
    B1 = 1000,
    B2 = 100
  )

  # Second derivative
  result_second <- compute_mean_deriv(
    Y = df[, Y],
    X = X,
    bws = bws,
    degree = 4,
    gradient.vec = 2,  # second derivative
    B1 = 1000,
    B2 = 100
  )

  # Append results (one row holds both the first and second derivative results)
  results <- rbind(results, data.frame(
    label = Y,
    First_Mean = result_first$mean_derivative,
    First_P = result_first$p_value,
    First_Sig = result_first$p_value < 0.05,
    Second_Mean = result_second$mean_derivative,
    Second_P = result_second$p_value,
    Second_Sig = result_second$p_value < 0.05
  ))
}

# Save the results to a CSV file
write.csv(results, result_path("derivative_results.csv"), row.names = FALSE)
cat("Results saved to derivative_results.csv\n")
