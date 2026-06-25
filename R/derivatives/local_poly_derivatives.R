# ---------------------------------------------------------------------------
# Exploratory local-polynomial (generalized local polynomial regression)
# derivative estimation and diagnostic plots for a single condition.
#
# Companion to crs_spline_derivatives.R: fits npglpreg with CV-selected
# bandwidths and inspects the fitted smooth together with its first and second
# derivatives. Here the smooth variable is the contrast-coded performance
# score (BIS_Contrast). Run from the repository root.
# ---------------------------------------------------------------------------

library("crs")
library("doParabar")
library("foreach")

source("R/config.R")
source("R/derivatives/boot_der.R")

df <- read.csv(data_path("contrast_face.csv"))
vars <- colnames(df)
Y_names <- vars[10:369]
X <- df[c("BIS_Contrast", "Age", "Gender")]
# convert the two covariates to factors
X$Age = ordered(X$Age)
X$Gender = factor(X$Gender)

for (Y in Y_names)
{
    f <- formula(paste(Y, "~ BIS_Contrast + ordered(Age) + factor(Gender)"))

    # automatic bandwidth selection, using least-squares CV.
    model_bw <- npglpreg(formula = f, data = df, bwtype = "auto", degree = 4, nmulti = 5, cv = "bandwidth")
    bws = model_bw$bws

    # mean first-derivative p-value via the residual bootstrap (boot_der.R)
    compute_mean_deriv(Y = df[, Y], X = X, bws = bws, degree = 4, gradient.vec = 1, B1 = 1000, B2 = 100)

    # Diagnostic plots of the fitted smooth and its first / second derivatives.
    # Note: refit and assign `model` before plotting; npglpreg may require
    # manually re-attaching model$x / model$y as a workaround for a crs
    # package bug (see plot.npglpreg in np.regression.glp.R for details).
    model <- npglpreg(tydat = df[, Y], txdat = X, bwtype = "fixed", bws = bws, degree = 4, gradient.vec = 1)
    plot(model, deriv = 0, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
    plot(model, deriv = 1, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
    plot(model, deriv = 2, ci = TRUE, mean = TRUE, plot.errors.boot.num = 50, plot.errors.type = "quantiles", plot.behavior = "data")
}
