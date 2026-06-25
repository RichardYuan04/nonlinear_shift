# ---------------------------------------------------------------------------
# Plot a fitted GAM activation-performance curve for a single parcel.
#
# plot_gam() fits Activation ~ s(BIS, k, m) + Gender + Age, predicts the
# smooth across the observed BIS range (covariates held at median / modal
# level) and overlays the fitted curve and 95% CI on the raw scatter. Run
# from repo root. See data/README.md for the expected input columns.
# ---------------------------------------------------------------------------

library(mgcv)
library(ggplot2)
library(dplyr)

source("R/config.R")

place = read.csv(data_path("WM_Place_2bk_cleaned.csv"))
face  = read.csv(data_path("WM_Face_2bk_cleaned.csv"))


# The function fits a GAM model to the data and generates a plot of the fitted values along with the confidence intervals.
# k is the number of knots in the smooth term.
# m is the penalty on the order of derivatives of the smooth term.
# The default value is 2, which means the second derivative is penalized.


plot_gam <- function(data, y_var, smooth_var, 
           covariates, knots = 3, m = 2,
           set_fx = FALSE, print_summary = FALSE, print_smooth_parameters = FALSE) {
  # Data preprocessing
  y <- data[, y_var]
  smooth_var_data <- data[, smooth_var]
  covariates_data <- data[, covariates]
  
  # Remove missing values
  non_na_index <- complete.cases(data[, c(y_var, smooth_var, covariates)])
  data <- data[non_na_index, ]
  
  # Construct model formula
  model_formula <- as.formula(sprintf("%s ~ s(%s, k = %s, m = %s, fx = %s) + %s", y_var, smooth_var, knots, m, set_fx, paste(covariates, collapse = " + ")))
  
  # Fit GAM model
  gam_model <- gam(model_formula, method = "REML", data = data)
  
  # Print anything?
  if (print_summary) {
  print(summary(gam_model))
  }
  
  if (print_smooth_parameters) {
  smooth_summary <- summary(gam_model)$s.table
  print("Smooth Parameters:")
  print(smooth_summary)
  }
  
  # Generate plot data
  np <- 1000 # number of predicted values
  df <- gam_model$model
  theseVars <- attr(gam_model$terms, "term.labels")
  varClasses <- attr(gam_model$terms, "dataClasses")
  thisResp <- as.character(gam_model$terms[[2]])
  
  # Line plot with no interaction
  thisPred <- data.frame(init = rep(0, np))
  
  for (v in seq_along(theseVars)) {
  thisVar <- theseVars[[v]]
  thisClass <- varClasses[thisVar]
  if (thisClass == "character") {
    df[, thisVar] <- as.factor(df[, thisVar])
    thisClass <- "factor"
  }
  
  if (thisVar == smooth_var) {
    thisPred[, smooth_var] <- seq(min(df[, smooth_var], na.rm = TRUE), max(df[, smooth_var], na.rm = TRUE), length.out = np)
  } else {
    tab.tmp <- table(df[, thisVar])
    levelact <- which.max(tab.tmp)
    switch(thisClass,
       "numeric" = { thisPred[, thisVar] <- median(df[, thisVar]) },
       "factor" = { thisPred[, thisVar] <- levels(df[, thisVar])[[levelact]] },
       "ordered" = { thisPred[, thisVar] <- levels(df[, thisVar])[[levelact]] }
    )
  }
  }
  pred <- thisPred %>% dplyr::select(-init)
  
  p <- data.frame(predict(gam_model, pred, se.fit = TRUE))
  pred <- cbind(pred, p)
  pred$selo <- pred$fit - 1.96 * pred$se.fit
  pred$sehi <- pred$fit + 1.96 * pred$se.fit
  pred$fit.C <- scale(pred$fit, center = TRUE, scale = FALSE) # zero-centered fitted values
  pred$fit.Z <- scale(pred$fit, center = TRUE, scale = TRUE) # z-scored fitted values
  pred$fit.floor <- pred$fit - pred$fit[1] # fitted values minus the initial fit
  pred$fit.ratio <- pred$fit / pred$fit[1] # fitted values divided by the initial fit
  pred$selo.C <- pred$fit - 1.96 * pred$se.fit - mean(pred$fit)
  pred$sehi.C <- pred$fit + 1.96 * pred$se.fit - mean(pred$fit)
  pred[, thisResp] <- 1
  
  # Plot the data
  ggplot(data, aes_string(x = smooth_var, y = y_var)) +
  geom_point(alpha = 0.2) +
  geom_line(data = pred, aes_string(x = smooth_var, y = "fit"), color = "red") +
  geom_ribbon(data = pred, aes_string(x = smooth_var, ymin = "selo", ymax = "sehi"), alpha = 0.2) +
  labs(x = smooth_var, y = y_var) +
  #ylim(-0.5, 0.5)+
  theme_minimal()
}