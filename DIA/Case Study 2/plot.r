# Configuration
#------------------------------------------------------------------------------
# Set analysis mode: "DIANN" or "Spectronaut"
ANALYSIS_MODE <- "Spectronaut"  # Change this to switch between DIANN and Spectronaut

# Set paths based on analysis mode
METRICS_DIR <- if(ANALYSIS_MODE == "DIANN") "metrics_DIANN" else "metrics_Spectronaut"
PLOT_DIR <- if(ANALYSIS_MODE == "DIANN") "Plot_files/DIANN" else "Plot_files/Spectronaut"

# Analysis parameters
FDR <- 0.05  # FDR threshold for plots

# Create output directory if it doesn't exist
if (!dir.exists(PLOT_DIR)) dir.create(PLOT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Load required libraries
suppressPackageStartupMessages({
    library(stringr)     
    library(dplyr)       
    library(ggplot2)     
    library(ggstatsplot) 
    library(ggsignif)   
})

# Load and clean up metric data
# Find the most recent metrics file
metrics_files <- list.files(METRICS_DIR, pattern = paste0("metrics_", ANALYSIS_MODE, "_FDR", FDR), full.names = TRUE)
latest_file <- metrics_files[which.max(file.info(metrics_files)$mtime)]

if (is.null(latest_file) || length(latest_file) == 0) {
    stop(sprintf("No metrics file found for %s analysis with FDR %.2f", ANALYSIS_MODE, FDR))
}

# Read the CSV file containing performance metrics
message(sprintf("Loading metrics from: %s", latest_file))
metrics <- read.csv(latest_file)

# Clean column names: remove ".rds" suffix
colnames(metrics) <- str_remove(colnames(metrics), ".rds")

# Further clean: remove "_LFQ" suffix from column names
colnames(metrics) <- str_remove(colnames(metrics), "_LFQ")

# Filter the metrics to keep only specific evaluation scores
metrics <- metrics[metrics$X %in% c(
  "F1_Score", 
  "nMCC", 
  "G_Mean", 
  "Fowlkes_Mallows_Index", 
  "Balanced_Accuracy", 
  "pauc"
), ]


# Reshape data to long format for plotting
# Convert data to long format: columns for variable (method+experiment) and value (score)
metrics.melt <- reshape2::melt(metrics)

# Extract DEA method (e.g., Limma, ROTS) from variable column by splitting at first underscore
metrics.melt$DEA <- str_split_fixed(metrics.melt$variable, fixed("_"), 2)[,1]

# Extract experiment name from variable column (removes everything after first underscore)
metrics.melt$exp <- str_split_fixed(metrics.melt$variable, fixed("_"), 2)[,2]
metrics.melt$exp <- sub("_.*", "", metrics.melt$exp)  # Clean any remaining suffix

# Loop through each metric and generate plots
for(score in unique(metrics.melt$X)){
  
  # Filter data for current metric
  df <- metrics.melt[metrics.melt$X == score, ]
  
  # Replace missing values with 0 to avoid plotting errors
  df[is.na(df)] <- 0
  
  # Create violin/boxplot-style comparison using ggbetweenstats
  p <- ggbetweenstats(
    data = df,
    x = DEA,                 # X-axis: DEA method
    y = value,               # Y-axis: metric score
    type = "robust",         # Use robust ANOVA (median-based)
    centrality.type = "nonparametric", # Use nonparametric statistics
    p.adjust.method = "BH",  # Adjust p-values using Benjamini-Hochberg
    pairwise.display = "none", # Don't display pairwise tests
    bf.message = FALSE,        # Disable Bayes Factor message
    results.subtitle = FALSE,  # Hide test result subtitle
    package = "ggsci",         # Use Lancet color palette
    palette = "lanonc_lancet",
    ylab = score,              # Y-axis label is the current score name
    centrality.label.args = list(
      size = 4.5,
      segment.linetype = 1,
      nudge_x = 0.25
    ),
    ggplot.component = list(  # Customize theme
      theme(
        axis.text.x = element_text(color = "black", face = "bold", size = 13),
        axis.text.y = element_text(color = "black", face = "bold", size = 13),
        axis.text = element_text(color = "black", face = "bold", size = 13),
        axis.title.y.left  = element_text(color = "black", face = "bold", size = 17),
        axis.title.x.bottom = element_text(color = "black", face = "bold", size = 17),
        text = element_text(color = "black", face = "bold"),
        legend.title = element_text(color = "black", face = "bold", size = 13.5),
        legend.text = element_text(color = "black", face = "bold", size = 13.5),
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank()
      )
    ),
    ggsignif.args = list(     # Appearance settings for significance bars (not shown since pairwise.display is none)
      textsize = 25,
      tip_length = 0.01,
      color = "black",
      vjust = 0.5,
      size = 0.2
    )
  ) + scale_x_discrete(labels = levels(as.factor(df$DEA)))  # Fix x-axis labels

  # Save plot to file with name indicating analysis mode, metric, and FDR level
  plot_filename <- file.path(PLOT_DIR, sprintf("%s_%s_FDR%.2f.jpg", score, ANALYSIS_MODE, FDR))
  
  ggsave(
    filename = plot_filename,
    plot = p,
    dpi = 300,
    width = 8.5,
    height = 5.5
  )
  
  message(sprintf("Saved plot: %s", plot_filename))
}

# Print completion message
message(sprintf("\nAnalysis completed for %s mode", ANALYSIS_MODE))
message(sprintf("Plots saved in: %s", normalizePath(PLOT_DIR)))
