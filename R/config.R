# ---------------------------------------------------------------------------
# Shared configuration for the R analysis scripts.
#
# All scripts are meant to be run from the repository root. Input CSVs (the
# subject-level, HCP-derived data described in data/README.md) are read from
# DATA_DIR; generated outputs are written to RESULTS_DIR. Both can be
# overridden with environment variables so no absolute paths are hard-coded.
# ---------------------------------------------------------------------------

DATA_DIR    <- Sys.getenv("NLS_DATA_DIR",    unset = "data")
RESULTS_DIR <- Sys.getenv("NLS_RESULTS_DIR", unset = "results")

if (!dir.exists(RESULTS_DIR)) {
  dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)
}

# Convenience helpers.
data_path   <- function(...) file.path(DATA_DIR, ...)
result_path <- function(...) file.path(RESULTS_DIR, ...)
