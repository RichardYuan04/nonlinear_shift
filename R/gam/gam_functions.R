# ---------------------------------------------------------------------------
# Core GAM helper functions (mgcv) shared across the analysis scripts.
#
# These adapt the modelling utilities from Sydnor et al. (2021) for the
# region-wise activation ~ s(BIS) + covariates models used in this project.
# Provides routines to: fit the performance smooth and extract its statistics
# (partial R-squared, ANOVA p-value), draw posterior smooths, test a covariate
# of interest, fit a factor-smooth interaction, and compute first/second
# smooth derivatives. source() this file before calling the functions.
# ---------------------------------------------------------------------------

library(readxl)
library(tidyr)
library(mgcv)
library(psych)
library(dplyr)
library(ecostats)
library(ggplot2)
library(gratia)
library(stringr)


#gam.fit.smooth: fit GAM with smooth variable and return the statistics of the smooth variable.
gam.fit.smooth <- function(dataset, region, smooth_var, covariates, knots = 3, m = 4, set_fx = FALSE, stats_only = FALSE) {
  
  # Fit the gam with smooth variable
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, m = %s, fx = %s) + %s", region, smooth_var, knots, m, set_fx, paste(covariates, collapse = " + ")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # GAM statistics
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  # Fit the gam without smooth variable (null model)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, paste(covariates, collapse = " + ")))
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  # Calculate partial R squared
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  
  # Anova test for smooth variable
  anova.smooth.pvalue <- anova.gam(gam.nullmodel, gam.model, test = 'Chisq')$`Pr(>Chi)`[2]
  
  # Prepare results
  stats.results <- cbind(region, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue)
  
  if (stats_only == TRUE) {
    return(stats.results)
  } else {
    return(stats.results)
  }
}


#gam.fit.smooth: fit GAM with smooth variable and return the statistics of the smooth variable.
gam.fit.smooth <- function(dataset, region, smooth_var, covariates, knots = 3, set_fx = FALSE, stats_only = FALSE) {
  
  # Fit the gam with smooth variable
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, paste(covariates, collapse = " + ")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # GAM statistics
  gam.smooth.F <- gam.results$s.table[3]
  gam.smooth.pvalue <- gam.results$s.table[4]
  
  #Get derivatives of the smooth function using finite differences
  derv <- derivatives(gam.model, data = gam.data, term = sprintf('s(%s)',smooth_var), 
                      interval = "simultaneous", unconditional = F)
  
  # Fit the gam without smooth variable (null model)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, paste(covariates, collapse = " + ")))
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  # Calculate partial R squared
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  
  mean.derivative <- mean(derv$.derivative)
  if(mean.derivative < 0){ #if the average derivative is less than 0, make the effect size estimate negative
    partialRsq <- partialRsq*-1}
  
  # Anova test for smooth variable
  anova.smooth.pvalue <- anova.gam(gam.nullmodel, gam.model, test = 'Chisq')$`Pr(>Chi)`[2]
  
  # Prepare results
  stats.results <- cbind(region, gam.smooth.F, gam.smooth.pvalue, partialRsq, anova.smooth.pvalue)
  
  if (stats_only == TRUE) {
    return(stats.results)
  } else {
    return(stats.results)
  }
}




#fit.gam: fit GAM with smooth variable and return the statistics of the smooth variable.
fit_gam <- function(data, region, smooth_var, covariates, knots = 3, m = 4, set_fx = FALSE, stats_only = TRUE) {
  y <- data[, region]
  smooth_var_data <- data[, smooth_var]
  covariates_data <- data[, covariates]
  
  non_na_index <- complete.cases(data[, c(region, smooth_var, covariates)])
  data <- data[non_na_index, ]
  

  model_formula <- as.formula(sprintf("%s ~ s(%s, k = %s, m = 4, fx = %s) + %s", region, smooth_var, knots, m, set_fx, paste(covariates, collapse = " + ")))
  

  gam.model <- gam(model_formula, method = "REML", data = data)
  

  gam_results <- summary(gam.model)
  gam_smooth_F <- gam_results$s.table[1, 3]
  gam_smooth_pvalue <- gam_results$s.table[1, 4]
  
  # Fit the gam without smooth variable (null model)
  nullmodel <- as.formula(sprintf("%s ~ %s", region, paste(covariates, collapse = " + ")))
  gam.nullmodel <- gam(nullmodel, method = "REML", data = data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  # Calculate partial R squared
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  
  ##Full versus reduced model anova p-value
  anova.smooth.pvalue <- anova.gam(gam.nullmodel,gam.model,test='Chisq')$`Pr(>Chi)`[2]
  
  # calculate the correlation between residuals
  # res <- residuals(gam.model)
  # PCorr_Test <- corr.test(res, data[,smooth_var], method = "pearson")
  # correstimate <- as.numeric(PCorr_Test$r)
  # corrp <- as.numeric(PCorr_Test$p)
  
  stats_results <- cbind(region, smooth_var, gam_smooth_F, gam_smooth_pvalue, 
                         partialRsq, anova.smooth.pvalue)
  data_results <- list()
  data_results[[1]] <- as.data.frame(stats_results)
  
  if (stats_only == TRUE) {
    return(stats_results)
  } else {
    return(data_results)
  }
}

#gam.posterior.smooths: fit GAM with smooth variable and return posterior smooths for the smooth variable.
gam.posterior.smooths <- function(dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_draws = FALSE) {
  
  # Set parameters
  npd <- as.numeric(draws)
  np <- as.numeric(increments)
  EPS <- 1e-07
  UNCONDITIONAL <- FALSE
  
  # Fit the gam
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, paste(covariates, collapse = " + ")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # Extract gam input data
  df <- gam.model$model
  
  # Create a prediction data frame
  thisPred <- data.frame(init = rep(0, np))
  
  theseVars <- attr(gam.model$terms, "term.labels")
  varClasses <- attr(gam.model$terms, "dataClasses")
  thisResp <- as.character(gam.model$terms[[2]])
  
  for (v in c(1:length(theseVars))) {
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisClass == "character") {
      df[, thisVar] <- as.factor(df[, thisVar])
      thisClass <- "factor"
    } #check if the column are str, then reformat it into factor
    if (thisVar == smooth_var) {
      thisPred[, smooth_var] <- seq(min(df[, smooth_var], na.rm = TRUE), max(df[, smooth_var], na.rm = TRUE), length.out = np)
    } else {
      switch(thisClass,
             "numeric" = { thisPred[, thisVar] <- median(df[, thisVar]) },
             "factor" = { thisPred[, thisVar] <- levels(df[, thisVar])[[1]] },
             "ordered" = { thisPred[, thisVar] <- levels(df[, thisVar])[[1]] }
      )
    }
  }
  pred <- thisPred %>% select(-init)
  
  # Estimate posterior smooth functions (fitted values) from simulated GAM posterior distribution
  Vb <- vcov(gam.model, unconditional = UNCONDITIONAL)
  sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb)
  X0 <- predict(gam.model, newdata = pred, type = "lpmatrix")
  predicted.smooth.values <- X0 %*% t(sims)
  colnames(predicted.smooth.values) <- sprintf("draw%s", seq(from = 1, to = npd))
  predicted.smooth.values <- cbind(as.numeric(pred[, smooth_var]), predicted.smooth.values)
  colnames(predicted.smooth.values)[1] <- sprintf("%s", smooth_var)
  predicted.smooth.values <- as.data.frame(predicted.smooth.values)
  
  # Smooth minimum/maximum values and credible intervals
  max.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is largest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x) {
                       round(predicted.smooth.values[, smooth_var][which.max(x)], 2)
                     }))
  max.y.range <- t(max.y.range)
  max.y <- median(max.y.range) #median value
  max.y.CI <- quantile(max.y.range, probs = c(0.025, 0.975))
  max.y.CI.lower <- max.y.CI[[1]]
  max.y.CI.upper <- max.y.CI[[2]]
  
  min.y.range <- predicted.smooth.values %>% #the value of smooth_var when y is lowest for each draw
    summarise(across(contains("draw"),
                     .fns = function(x) {
                       round(predicted.smooth.values[, smooth_var][which.min(x)], 2)
                     }))
  min.y.range <- t(min.y.range)
  min.y <- median(min.y.range)
  min.y.CI <- quantile(min.y.range, probs = c(0.025, 0.975))
  min.y.CI.lower <- min.y.CI[[1]]
  min.y.CI.upper <- min.y.CI[[2]]
  
  if (return_draws == TRUE) {
    return(predicted.smooth.values)
  } else {
    smooth.features <- list(region, max.y, max.y.CI.lower, max.y.CI.upper, min.y, min.y.CI.lower, min.y.CI.upper)
    names(smooth.features) <- c("region", sprintf("%s at max y", smooth_var), "max y credible interval lower", "max y credible interval upper", sprintf("%s at min y", smooth_var), "min y credible interval lower", "min y credible interval upper")
    return(smooth.features)
  }
}





#gam.fit.covariate: fit GAM with smooth variable and covariate of interest, and return the covariate statistics.
gam.fit.covariate <- function(dataset, region, smooth_var, covariate.interest, covariates.noninterest, knots = 3, set_fx = FALSE) {
  
  # Fit the gam
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s + %s", region, smooth_var, knots, set_fx, covariate.interest, paste(covariates.noninterest, collapse = " + ")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # GAM statistics
  gam.cov.tvalue <- gam.results$p.table[2, 3]
  gam.cov.pvalue <- gam.results$p.table[2, 4]
  
  nullmodel <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, paste(covariates.noninterest, collapse = " + ")))
  gam.nullmodel <- gam(nullmodel, method = "REML", data = gam.data)
  gam.nullmodel.results <- summary(gam.nullmodel)
  
  anova.cov.pvalue <- anova.gam(gam.nullmodel, gam.model, test = 'Chisq')$`Pr(>Chi)`[2]
  if (is.na(anova.cov.pvalue)) {
    anova.cov.pvalue <- 1
  }
  
  sse.model <- sum((gam.model$y - gam.model$fitted.values)^2)
  sse.nullmodel <- sum((gam.nullmodel$y - gam.nullmodel$fitted.values)^2)
  partialRsq <- (sse.nullmodel - sse.model) / sse.nullmodel
  if (gam.cov.tvalue < 0) {
    partialRsq <- partialRsq * -1
  }
  
  results <- cbind(region, gam.cov.tvalue, gam.cov.pvalue, anova.cov.pvalue, partialRsq)
  return(results)
}





#gam.factorsmooth.interaction: fit GAM with smooth variable and factor interaction, and return the interaction statistics.
gam.factorsmooth.interaction <- function(dataset, region, smooth_var, int_var, covariates, knots, set_fx = FALSE) {
  
  # Convert character columns to factor
  dataset <- as.data.frame(lapply(dataset, function(col) {
    if (is.character(col)) {
      return(as.factor(col))
    } else {
      return(col)
    }
  }))
  
  # Fit the gam
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%1$s ~ s(%2$s, k = %3$s, fx = %4$s) + s(%2$s, by = %5$s, k = %3$s, fx = %4$s) + %6$s", region, smooth_var, knots, set_fx, int_var, paste(covariates, collapse = "+")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # GAM statistics
  gam.int.F <- gam.results$s.table[2, 3]
  gam.int.pvalue <- gam.results$s.table[2, 4]
  
  interaction.stats <- cbind(region, gam.int.F, gam.int.pvalue)
  return(interaction.stats)
}


#gam.derivatives: compute smooth derivatives for the main GAM model and also for different draws from the simulated posterior distribution
gam.derivatives <- function(dataset, region, smooth_var, covariates, knots, set_fx = FALSE, draws, increments, return_posterior_derivatives = TRUE){
  
  # Convert character columns to factor
  dataset <- as.data.frame(lapply(dataset, function(col) {
    if (is.character(col)) {
      return(as.factor(col))
    } else {
      return(col)
    }
  }))
  
  #Set parameters
  npd <- as.numeric(draws) #number of draws from the posterior distribution; number of posterior derivative sets estimated
  np <- as.numeric(increments) #number of smooth_var increments to get derivatives at
  EPS <- 1e-07 #finite differences
  UNCONDITIONAL <- FALSE #should we account for uncertainty when estimating smoothness parameters?
  
  #Fit the gam
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, paste(covariates, collapse = "+")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  #Extract gam input data
  df <- gam.model$model #extract the data used to build the gam, i.e., a df of y + predictor values 
  
  #Create a prediction data frame, used to estimate (posterior) model coefficients
  thisPred <- data.frame(init = rep(0,np)) 
  
  theseVars <- attr(gam.model$terms,"term.labels") #gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms,"dataClasses") #classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) #the measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { #fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[,smooth_var] = seq(min(df[,smooth_var],na.rm = T),max(df[,smooth_var],na.rm = T), length.out = np) #generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[,thisVar] = median(df[,thisVar])}, #make predictions based on median value
              "factor" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]}, #make predictions based on first level of factor 
              "ordered" = {thisPred[,thisVar] = levels(df[,thisVar])[[1]]} #make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init) #prediction df
  pred2 <- pred #second prediction df
  pred2[,smooth_var] <- pred[,smooth_var] + EPS #finite differences
  
  #Estimate smooth derivatives
  derivs <- derivatives(gam.model, select = sprintf('s(%s)',smooth_var),#derivative at 200 indices of smooth_var with a simultaneous CI
                        interval = "simultaneous", 
                        unconditional = UNCONDITIONAL, 
                        data = pred) #used to be 'newdata' argument
  derivs.fulldf <- derivs %>% dplyr::select(BIS, .derivative, .se, .lower_ci, .upper_ci)
  # derivs.fulldf <- derivs %>% 
  #   mutate(smooth_var = data[[smooth_var]]) %>%
  #   select(smooth_var, derivative, se, lower, upper)
  derivs.fulldf <- derivs.fulldf %>% mutate(significant = !(0 > .lower_ci & 0 < .upper_ci))
  derivs.fulldf$significant.derivative = derivs.fulldf$.derivative*derivs.fulldf$significant
  colnames(derivs.fulldf) <- c(sprintf("%s", smooth_var), "derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  #Estimate posterior smooth derivatives from simulated GAM posterior distribution
  if(return_posterior_derivatives == TRUE){
    Vb <- vcov(gam.model, unconditional = UNCONDITIONAL) #variance-covariance matrix for all the fitted model parameters (intercept, covariates, and splines)
    
    #use matrix transform when npd equals 1 to escape package setting error
    if(npd == 1){
      sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb) #simulate model parameters (coefficents) from the posterior distribution of the smooth based on actual model coefficients and covariance
    sims <- t(as.matrix(sims))}else{
      sims <- MASS::mvrnorm(npd, mu = coef(gam.model), Sigma = Vb)
    }
    
    # 
    
    X0 <- predict(gam.model, newdata = pred, type = "lpmatrix") #get matrix of linear predictors for pred
    X1 <- predict(gam.model, newdata = pred2, type = "lpmatrix") #get matrix of linear predictors for pred2
    Xp <- (X1 - X0) / EPS 
    posterior.derivs <- Xp %*% t(sims) #Xp * simulated model coefficients = simulated derivatives. Each column of posterior.derivs contains derivatives for a different draw from the simulated posterior distribution
    posterior.derivs <- as.data.frame(posterior.derivs)
    colnames(posterior.derivs) <- sprintf("draw%s",seq(from = 1, to = npd)) #label the draws
    posterior.derivs <- cbind(as.numeric(pred[,smooth_var]), posterior.derivs) #add smooth_var increments from pred df to first column
    colnames(posterior.derivs)[1] <- sprintf("%s", smooth_var) #label the smooth_var column
    posterior.derivs <- cbind(as.character(region), posterior.derivs) #add parcel label to first column
    colnames(posterior.derivs)[1] <- "label" #label the column
    posterior.derivs.long <- posterior.derivs %>% pivot_longer(contains("draw"), names_to = "draw",values_to = "posterior.derivative")
  } #np*npd rows, 3 columns (smooth_var, draw, posterior.derivative)
  
  if(return_posterior_derivatives == FALSE)
    return(derivs.fulldf)
  if(return_posterior_derivatives == TRUE)
    return(posterior.derivs.long)
}





#gam.second.derivative: calculate the second derivative of the smooth term in the GAM model.
gam.second.derivative <- function(dataset, region, smooth_var, covariates, knots, set_fx = FALSE, increments){
  
  # Convert character columns to factor
  dataset <- as.data.frame(lapply(dataset, function(col) {
    if (is.character(col)) {
      return(as.factor(col))
    } else {
      return(col)
    }
  }))
  
  # Set parameters
  np <- as.numeric(increments) # Number of smooth_var increments to get derivatives at
  EPS <- 1e-07 # Finite differences
  UNCONDITIONAL <- FALSE # Should we account for uncertainty when estimating smoothness parameters?
  
  # Fit the gam
  gam.data <- dataset
  modelformula <- as.formula(sprintf("%s ~ s(%s, k = %s, fx = %s) + %s", region, smooth_var, knots, set_fx, paste(covariates, collapse = "+")))
  gam.model <- gam(modelformula, method = "REML", data = gam.data)
  gam.results <- summary(gam.model)
  
  # Extract gam input data
  df <- gam.model$model # Extract the data used to build the gam, i.e., a df of y + predictor values 
  
  # Create a prediction data frame, used to estimate (posterior) model coefficients
  thisPred <- data.frame(init = rep(0, np)) 
  
  theseVars <- attr(gam.model$terms, "term.labels") # Gam model predictors (smooth_var + covariates)
  varClasses <- attr(gam.model$terms, "dataClasses") # Classes of the model predictors and y measure
  thisResp <- as.character(gam.model$terms[[2]]) # The measure to fit for each posterior draw
  for (v in c(1:length(theseVars))) { # Fill the prediction df with data 
    thisVar <- theseVars[[v]]
    thisClass <- varClasses[thisVar]
    if (thisVar == smooth_var) { 
      thisPred[, smooth_var] = seq(min(df[, smooth_var], na.rm = TRUE), max(df[, smooth_var], na.rm = TRUE), length.out = np) # Generate a range of np data points, from minimum of smooth term to maximum of smooth term
    } else {
      switch (thisClass,
              "numeric" = {thisPred[, thisVar] = median(df[, thisVar])}, # Make predictions based on median value
              "factor" = {thisPred[, thisVar] = levels(df[, thisVar])[[1]]}, # Make predictions based on first level of factor 
              "ordered" = {thisPred[, thisVar] = levels(df[, thisVar])[[1]]} # Make predictions based on first level of ordinal variable
      )
    }
  }
  pred <- thisPred %>% dplyr::select(-init) # Prediction df
  pred2 <- pred # Second prediction df
  pred2[, smooth_var] <- pred[, smooth_var] + EPS # Finite differences
  
  # Estimate second derivatives
  second_derivs <- derivatives(gam.model, select = sprintf('s(%s)', smooth_var), # Second derivative at 200 indices of smooth_var with a simultaneous CI
                               order = 2, 
                               interval = "simultaneous", 
                               unconditional = UNCONDITIONAL, 
                               data = pred) # Used to be 'newdata' argument
  second_derivs.fulldf <- second_derivs %>% dplyr::select(BIS, .derivative, .se, .lower_ci, .upper_ci)
  second_derivs.fulldf <- second_derivs.fulldf %>% mutate(significant = !(0 > .lower_ci & 0 < .upper_ci))
  second_derivs.fulldf$significant.derivative = second_derivs.fulldf$.derivative * second_derivs.fulldf$significant
  colnames(second_derivs.fulldf) <- c(sprintf("%s", smooth_var), "second_derivative", "se", "lower", "upper", "significant", "significant.derivative")
  
  return(second_derivs.fulldf)
}

