# --------------------------------------------------------------------------------------------------------------------------------------------
# This implementation is based on:
# Racine, J. (1997). Consistent Significance Testing for Nonparametric Regression, Journal of Business & Economic Statistics, 15:3, 369-378
# and the compute.bootstrap.errors() function in the crs package (in /R/np.regression.glp.R): 
# Racine J., Nie Z. (2024). crs: Categorical Regression Splines. R package version 0.15-38, <https://github.com/JeffreyRacine/R-Package-crs>.
# --------------------------------------------------------------------------------------------------------------------------------------------

compute_mean_deriv <- function(
        Y,             # numeric vector of responses (length n)
        X,             # data.frame of covariates
        bws,           # bandwidth(s) for npglpreg
        degree,        # polynomial degree for npglpreg
        gradient.vec,  # order of derivative(s) to extract (integer or vector)
        B1 = 1000,      # number of outer bootstrap repetitions
        B2 = 100,       # number of inner bootstrap repetitions
        n_cores = NULL
) {
    n <- length(Y)
    
    # ----------------------------------------------------------------------------
    # 1) Fit the original model and compute observed mean derivative
    # ----------------------------------------------------------------------------
    cat(sprintf("1. Fit original model... "))
    
    orig.fit <- npglpreg(
        tydat = Y,
        txdat = X,
        bwtype = "fixed",
        bws = bws,
        degree = degree,
        gradient.vec = gradient.vec
    )
    lambda_obs <- mean(orig.fit$gradient)
    
    cat(sprintf("Done.\n"))
    
    # ----------------------------------------------------------------------------
    # 2) Inner bootstrap: estimate SE of the mean derivative via resampling pairs
    # ----------------------------------------------------------------------------
    estimate_se <- function(y_vec, x_mat) {
        boot_means <- numeric(B2)
        
        for (b in seq_len(B2)) {
            idx <- sample.int(n, n, replace = TRUE)
            y_b <- y_vec[idx]
            x_b <- x_mat[idx, , drop = FALSE]
            fit_b <- npglpreg(
                tydat = y_b,
                txdat = x_b,
                bwtype = "fixed",
                bws = bws,
                degree = degree,
                gradient.vec = gradient.vec
            )
            boot_means[b] <- mean(fit_b$gradient)
        }
        
        # standard deviation of the bootstrap means
        sd(boot_means)
    }
    
    cat(sprintf("2. Estimate SE of the mean derivative via bootstrap (B2 = %d)... ", B2))
    
    SE_obs <- estimate_se(Y, X)
    t_obs <- lambda_obs / SE_obs  # pivot
    
    cat(sprintf("Done.\n"))
    
    # ----------------------------------------------------------------------------
    # 3) Outer residual bootstrap under H0: mean derivative = 0
    # ----------------------------------------------------------------------------
    
    cat(sprintf("3. Outer residual bootstrap (B1 = %d)...\n", B1))
    
    # 3a) Extract and center residuals under the null fit
    X_null <- X
    X_null[, "BIS_Contrast"] <- mean(X_null[, "BIS_Contrast"])
    Y_hat = predict(orig.fit, new_data = X_null)
    #Y_hat = orig.fit$fitted.values
    residuals0 <- Y - Y_hat
    residuals0 <- residuals0 - mean(residuals0)
    
    # 3b) Set up parallel backend
    if (is.null(n_cores)) {
        n_cores <- parallel::detectCores()
    }
    backend <- parabar::start_backend(cores = n_cores - 1, cluster_type = "fork", backend_type = "async") # change cluster_type to "psock" if on Windows
    registerDoParabar(backend)
    
    # 3c) Parallel calculation of pivoted statistics
    t_star <- foreach(i = seq_len(B1), .export = c("Y_hat", "residuals0", "X", "bws", "degree", "gradient.vec", "estimate_se", "n"), .combine = 'c', .packages = 'crs') %dopar% {
        # generate bootstrap response
        eps_star <- sample(residuals0, n, replace = TRUE)
        Y_star <- Y_hat + eps_star
        
        # refit and compute mean derivative
        fit_i <- npglpreg(
            tydat = Y_star,
            txdat = X,
            bwtype = "fixed",
            bws = bws,
            degree = degree,
            gradient.vec = gradient.vec
        )
        lambda_i <- mean(fit_i$gradient)
        
        # estimate SE for this replicate
        SE_i <- estimate_se(Y_star, X)
        
        # pivot
        lambda_i / SE_i
    }
    stop_backend(backend)
    
    cat(sprintf("Done.\n"))
    
    # ----------------------------------------------------------------------------
    # 4) Compute two-sided p-value for H0: mean derivative = 0
    # ----------------------------------------------------------------------------
    p_value <- mean(abs(t_star) >= abs(t_obs))
    
    # ----------------------------------------------------------------------------
    # Return results
    # ----------------------------------------------------------------------------
    list(
        mean_derivative = lambda_obs,
        SE_mean_derivative = SE_obs,
        t_observed = t_obs,
        t_star = t_star,
        p_value = p_value
    )
}
