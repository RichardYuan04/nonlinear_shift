# ---------------------------------------------------------------------------
# Cortical surface maps of parcel response type (ggseg + Glasser HCP-MMP).
#
# Each parcel is classified into one of six response patterns from the sign
# and magnitude of its mean first and second derivatives (linear / concave /
# convex x increase / decrease), then painted onto the cortical surface.
# Input CSVs are the per-parcel derivative results produced by
# bootstrap_derivatives.R (one file per condition). Figures are written to
# RESULTS_DIR. Run from the repository root.
# ---------------------------------------------------------------------------

library(dplyr)
library(ggseg)
library(ggsegGlasser)
library(ggplot2)
library(stringr)

source("R/config.R")

# Load the Glasser atlas labels
all_brain_areas <- glasser$data$label #get parcels' labels

# Per-condition derivative result files (see bootstrap_derivatives.R output).
data_paths <- list(
  WM_Body_2bk  = data_path("WM_Body_2bk_Sig_gam_der.csv"),
  WM_Place_2bk = data_path("WM_Place_2bk_Sig_gam_der.csv"),
  WM_Face_2bk  = data_path("WM_Face_2bk_Sig_gam_der.csv"),
  WM_Tool_2bk  = data_path("WM_Tool_2bk_Sig_gam_der.csv")
)


# This loop generates a brain map for each dataset, coloring regions based on their type.
# It uses the ggseg and ggsegGlasser packages to visualize the brain regions and their types.

for (name in names(data_paths)) {
  # Read the data frame, add error handling
  df <- tryCatch({
    read.csv(data_paths[[name]])
  }, error = function(e) {
    warning(paste("Failed to read file:", data_paths[[name]]))
    return(NULL)
  })
  
  if (is.null(df)) next  # If reading fails, skip the current loop
  
  # Add rows for missing brain regions, fill with NA or default values
  df_full <- data.frame(label = all_brain_areas) %>%
    left_join(df, by = "label") %>%
    mutate(
      First_Mean = ifelse(is.na(First_Mean), 0, First_Mean),
      Second_Mean = ifelse(is.na(Second_Mean), 0, Second_Mean)
    )
  
  # Assign region types based on the values of the first and second derivatives
  df_full <- df_full %>%
    mutate(
      region_type = case_when(
        First_Mean > 0 & Second_Mean > 0.05 ~ "Convex Increase",
        First_Mean > 0 & Second_Mean < -0.05 ~ "Concave Increase",
        First_Mean < 0 & Second_Mean > 0.05 ~ "Convex Decrease",
        First_Mean < 0 & Second_Mean < -0.05 ~ "Concave Decrease",
        Second_Mean >= -0.05 & Second_Mean <= 0.05 & First_Mean > 0 ~ "Linear Increase",
        Second_Mean >= -0.05 & Second_Mean <= 0.05 & First_Mean < 0 ~ "Linear Decrease",
        TRUE ~ "No Response"
      )
    )
  
  keyword <- str_extract(name, "(?<=_)[A-Za-z]+(?=_)")  # Extract the part between two underscores
  
  # Define color mapping
  region_colors <- c(
    "Convex Increase" = "red",
    "Concave Increase" = "orange",
    "Convex Decrease" = "green",
    "Concave Decrease" = "blue",
    "Linear Increase" = "pink",
    "Linear Decrease" = "cyan",
     "No Response" = "gray30"
  )
  
  # Create brain map
  p <- ggseg(
    df_full,
    atlas = glasser,
    mapping = aes(
      fill = region_type,  # Fill color based on region type
      colour = I("gray")   # Set border color to gray uniformly
    ),
    position = "stacked",
    size = 0.3  # Set border size uniformly
  ) +
    theme_void() +
    scale_fill_manual(values = region_colors) +  # Use custom color mapping
    labs(
      fill = "Region Type"
    ) +
    ggtitle(keyword) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.text = element_text(family = "Arial", color = c("black"))
    )
  
  # Save the image
  ggsave(
    filename = result_path(paste0(name, "_region_type_brainmap.jpg")),
    plot = p,
    device = "jpg",
    dpi = 500,
    width = 8,
    height = 6
  )
}




# Sample single brain plot for different region types


library(ggseg)
library(ggsegGlasser)
library(dplyr)
library(stringr)

# Data path for sig_gam_der:
data_paths <- list(
  WM_Body_2bk  = data_path("WM_Body_2bk_Sig_gam_der.csv"),
  WM_Place_2bk = data_path("WM_Place_2bk_Sig_gam_der.csv"),
  WM_Face_2bk  = data_path("WM_Face_2bk_Sig_gam_der.csv"),
  WM_Tool_2bk  = data_path("WM_Tool_2bk_Sig_gam_der.csv")
)

# Define function: Generate brain map for a single brain region

# Define target brain region (e.g., left hemisphere V1)
target_label <- "rh_R_a9-46v"  # Replace with your target brain region label


# This loop generates a brain map for a single brain region with specified datasets and region type.
# It uses the ggseg and ggsegGlasser packages to visualize the brain regions and their types.
# Need to manually modify the border color in the function!

for (name in names(data_paths)) {
  # Read the data frame, add error handling
  df <- tryCatch({
    read.csv(data_paths[[name]])
  }, error = function(e) {
    warning(paste("Failed to read file:", data_paths[[name]]))
    return(NULL)
  })
  
  if (is.null(df)) next  # If reading fails, skip the current loop
  
  # Add rows for missing brain regions, fill with NA or default values
  df_full <- data.frame(label = all_brain_areas) %>%
    left_join(df, by = "label") %>%
    mutate(
      First_Mean = ifelse(is.na(First_Mean), 0, First_Mean),
      Second_Mean = ifelse(is.na(Second_Mean), 0, Second_Mean)
    )
  
  # Assign region types based on the values of the first and second derivatives
  df_full <- df_full %>%
    mutate(
      region_type = case_when(
        First_Mean > 0 & Second_Mean > 0.05 ~ "Convex Increase",
        First_Mean > 0 & Second_Mean < -0.05 ~ "Concave Increase",
        First_Mean < 0 & Second_Mean > 0.05 ~ "Convex Decrease",
        First_Mean < 0 & Second_Mean < -0.05 ~ "Concave Decrease",
        Second_Mean >= -0.05 & Second_Mean <= 0.05 & First_Mean > 0 ~ "Linear Increase",
        Second_Mean >= -0.05 & Second_Mean <= 0.05 & First_Mean < 0 ~ "Linear Decrease",
        TRUE ~ "No Response"
      )
    ) %>%
    # Override colors: Only keep the color for the target brain region, others set to white
    mutate(
      region_type = ifelse(label == target_label, region_type, "Other")
    )
  
  # Define color mapping (add "Other" corresponding to white)
  region_colors <- c(
    "Convex Increase" = "red",
    "Concave Increase" = "orange",
    "Convex Decrease" = "green",
    "Concave Decrease" = "blue",
    "Linear Increase" = "pink",
    "Linear Decrease" = "cyan",
    "No Response" = "gray30",
    "Other" = "white"  # Other brain regions set to white
  )
  
  # Extract metadata for the target brain region (hemisphere and side view)
  target_meta <- glasser$data %>% filter(label == target_label)
  target_hemi <- unique(target_meta$hemi)   # e.g., "left"
  target_side <- unique(target_meta$side)   # e.g., "lateral"
  
  # Create brain map
  p <- ggseg(
    df_full,
    atlas = glasser,
    mapping = aes(
      fill = region_type,  # Fill color based on region type
      colour = I("blue")   # Set border color to blue uniformly
    ),
    hemisphere = target_hemi,  # Explicitly specify hemisphere
    side = target_side,        # Explicitly specify side view (lateral/medial)
    position = "stacked",
    size = 0.3
  ) +
    theme_void() +
    scale_fill_manual(values = region_colors) +  # Use custom color mapping
    labs(
      fill = "Region Type"
    ) +
    ggtitle(paste(target_label, "-", name)) +  # Title includes target brain region and dataset name
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      legend.text = element_text(family = "Arial", color = c("black"))
    )
  
  # Save the image
  ggsave(
    filename = result_path(paste0(target_label, "_", name, "_brainmap.jpg")),
    plot = p,
    device = "jpg",
    dpi = 500,
    width = 8,
    height = 6
  )
}




# Example usage
WM_Place_2bk <- read.csv(data_path("WM_Place_2bk_Sig_gam_der.csv"))
plot_single_region(input_label = "lh_L_V1", dataset_name = "WM_Place_2bk")
