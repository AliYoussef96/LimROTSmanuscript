# Load required libraries

library(stringr)     
library(dplyr)       
library(ggplot2)     
library(ggstatsplot) 
library(ggsignif)   

# Load and clean up metric data


# Read the CSV file containing various performance metrics for different DEA methods and experiments
metrics <- read.csv("metrics_Maxquant//FDR_0.05")

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

# Set FDR value to be used in file names
FDR <- 0.05
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
    pairwise.display = "none", # Don’t display pairwise tests
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
        axis.text.x = element_text(color = "black", face = "bold", size = 12),
        axis.text.y = element_text(color = "black", face = "bold", size = 12),
        axis.text = element_text(color = "black", face = "bold", size = 12),
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

  # Save plot to file with name indicating the metric and FDR level
  ggsave(
    filename = paste0("metrics_Maxquant//", score, "_FDR", FDR, ".jpg"),
    plot = p,
    dpi = 300,
    width = 8.5,
    height = 5.5
  )
}
