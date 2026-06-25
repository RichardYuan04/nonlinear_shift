# ---------------------------------------------------------------------------
# Bootstrap significance testing of the mean first and second derivatives of
# the GAM performance smooth, per parcel (mgcv + gratia).
#
# For each parcel we fit Activation ~ s(BIS, k = 3) + Gender + Age, then
# bootstrap a pivoted statistic (mean derivative / its bootstrap SE) to obtain
# a two-sided p-value for H0: mean derivative = 0. The sign of the mean first
# derivative gives the direction (increase / decrease); the sign of the mean
# second derivative gives the curvature (convex / concave). Output is one row
# per parcel written to RESULTS_DIR. Run from the repository root.
# ---------------------------------------------------------------------------

library(dplyr)
library(boot)
library(mgcv)
library(gratia)
library(tools)

source("R/config.R")

# --- Helper Functions ---

# Bootstrap for the second derivative
bootstrap_second_derivative <- function(data, indices, model, smooth_var, population_mean = 0, get_mean = FALSE) {
    boot_sample <- data[indices, ]
    derivative <- derivatives(model, order = 2, data = boot_sample, term = sprintf("s(%s)", smooth_var), type = "central")$.derivative
    derivative <- derivative - population_mean # demean if necessary
    if (get_mean) {
        return(list(mean(derivative), mean(derivative) * sqrt(length(derivative)) / sd(derivative)))
    }
    return(mean(derivative) * sqrt(length(derivative)) / sd(derivative))
}

# Bootstrap for the first derivative
bootstrap_first_derivative <- function(data, indices, model, smooth_var, population_mean = 0, get_mean = FALSE) {
    boot_sample <- data[indices, ]
    derivative <- derivatives(model, order = 1, data = boot_sample, term = sprintf("s(%s)", smooth_var), type = "central")$.derivative
    derivative <- derivative - population_mean # demean if necessary
    if (get_mean) {
        return(list(mean(derivative), mean(derivative) * sqrt(length(derivative)) / sd(derivative)))
    }
    return(mean(derivative) * sqrt(length(derivative)) / sd(derivative))
}

# Perform bootstrap resampling and calculate p-values
bootstrap_and_test <- function(df, model, dependent_var, smooth_var, nboot, bootstrap_func) {
    tryCatch({
        # Observed t-statistics and means
        obs_result <- bootstrap_func(df, 1:nrow(df), model, smooth_var, get_mean = TRUE)
        
        # Extract the mean values and observed t-statistics
        mean_all <- obs_result[[1]]
        obs_all <- obs_result[[2]]
        
        # Perform bootstrap resampling
        boot_all <- boot(df, statistic = function(data, indices) bootstrap_func(data, indices, model, smooth_var, population_mean = mean_all), R = nboot)
        
        # Compute p-values
        p_all <- (1 + sum(abs(boot_all$t) >= abs(obs_all))) / (nboot + 1)
        
        # Return results
        return(data.frame(
            Mean = mean_all,
            P = p_all,
            Sig = p_all < 0.05
        ))
    }, error = function(e) {
        cat("Error in ROI:", dependent_var, "\n", e$message, "\n")
        return(NULL)
    })
}

# --- Core Analysis Functions ---

# Analyze a single ROI for first-order and second-order derivatives
analyze_roi <- function(df, dependent_var, smooth_var, nboot) {
    cat(sprintf("Analyzing ROI: %s\n", dependent_var))
    
    # Fit the GAM model on the full dataset
    model <- gam(as.formula(paste(dependent_var, "~ s(", smooth_var, ", k = 3) + Gender + Age")), data = df, method = "REML")
    
    # First-order derivative
    first_order_results <- bootstrap_and_test(df, model, dependent_var, smooth_var, nboot, bootstrap_first_derivative)
    if (is.null(first_order_results)) {
        first_order_results <- data.frame(
            Mean = NA,
            P = NA,
            Sig = NA
        )
    }
    
    # Second-order derivative
    second_order_results <- bootstrap_and_test(df, model, dependent_var, smooth_var, nboot, bootstrap_second_derivative)
    if (is.null(second_order_results)) {
        second_order_results <- data.frame(
            Mean = NA,
            P = NA,
            Sig = NA
        )
    }
    
    # Combine results into a single row
    results <- data.frame(
        ROI = dependent_var,
        
        # First-order Results
        First_Mean = first_order_results$Mean,
        First_P = first_order_results$P,
        First_Sig = first_order_results$Sig,
        
        # Second-order Results
        Second_Mean = second_order_results$Mean,
        Second_P = second_order_results$P,
        Second_Sig = second_order_results$Sig,
        
        stringsAsFactors = FALSE
    )
    
    return(results)
}

# Run the analysis for a given file (file_name is resolved against DATA_DIR)
run_analysis <- function(file_name, smooth_var, nboot, seed) {
    set.seed(seed)
    file_path <- data_path(file_name)
    cat("Processing file:", file_path, "\n")

    # Load data
    df <- read.csv(file_path)
    
    # Extract ROI variables
    roi_vars <- colnames(df)[10:ncol(df)]
    results <- data.frame()
    
    # Analyze each ROI
    for (roi in roi_vars) {
        roi_result <- analyze_roi(df, roi, smooth_var, nboot)
        if (!is.null(roi_result)) results <- rbind(results, roi_result)
    }
    return(results)
}

# --- Main Script Execution ---
main <- function() {
    # CSV file names located in DATA_DIR (see data/README.md). Add the other
    # conditions (Face / Body / Tool) here to analyse them as well.
    file_names <- c(
        # "WM_Face_2bk_cleaned.csv",
        # "WM_Body_2bk_cleaned.csv",
        # "WM_Tool_2bk_cleaned.csv",
        "WM_Place_2bk_cleaned.csv"
    )
    smooth_var <- "BIS"
    nboot <- 2000
    seed <- 100

    for (file_name in file_names) {
        results <- run_analysis(file_name, smooth_var, nboot, seed)
        output_file <- result_path(paste0(file_path_sans_ext(basename(file_name)), "_analysis_results.csv"))
        write.csv(results, output_file, row.names = FALSE)
        cat("Results saved to:", output_file, "\n")
    }
}

# Run the main function
main()
