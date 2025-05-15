# Load necessary libraries
library(stringr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)

# Load and clean metric data from Spectronaut results
metrics <- read.csv("metrics_Spectronaut//FDR_0.05")
colnames(metrics) <- str_remove(colnames(metrics), ".rds")  # Remove ".rds" suffix from column names

# Filter to retain only selected evaluation metrics
metrics <- metrics[metrics$X %in% c(
  "F1_Score", "nMCC", "G_Mean", "Fowlkes_Mallows_Index",
  "Balanced_Accuracy", "pauc"
), ]

# Melt the dataframe into long format for easier plotting
metrics.melt <- reshape2::melt(metrics)

# Extract DEA method and experiment name from variable names
metrics.melt$DEA <- str_split_fixed(metrics.melt$variable, fixed("_"), 2)[, 1]
metrics.melt$exp <- str_split_fixed(metrics.melt$variable, fixed("_"), 2)[, 2]
metrics.melt$exp <- sub("_.*", "", metrics.melt$exp)  # Remove any extra suffix in `exp`

# Compute the median metric score for each combination of score type, DEA method, and experiment
median.DEA <- metrics.melt %>%
  group_by(X, DEA, exp) %>%
  summarise(median.DEA = median(value, na.rm = TRUE), .groups = "drop") %>%
  data.frame()

# Barplot: median values of each score per DEA method and experiment
ggplot(median.DEA, aes(x = median.DEA, y = X, fill = DEA)) + 
  geom_bar(stat = "identity", position = position_dodge()) + 
  theme_bw() + 
  facet_grid(. ~ exp)

# Set FDR threshold for labeling saved files
FDR <- 0.05

# Generate violin + significance plots per score type using ggbetweenstats
for (score in unique(metrics.melt$X)) {

  df <- metrics.melt[metrics.melt$X == score, ]
  df[is.na(df)] <- 0  # Replace NAs with 0

  p <- ggbetweenstats(
    data = df,
    x = DEA,
    y = value,
    type = "nonparametric",
    centrality.type = "nonparametric",
    p.adjust.method = "BH",
    pairwise.display = "s",
    bf.message = FALSE,
    results.subtitle = FALSE,
    package = "ggsci",
    palette = "lanonc_lancet",
    ylab = score,
    ggplot.component = list(
      theme(
        axis.text.x = element_text(color = "black", face = "bold", size = 7),
        axis.text.y = element_text(color = "black", face = "bold", size = 9),
        axis.text = element_text(color = "black", face = "bold", size = 7),
        text = element_text(color = "black", face = "bold"),
        axis.title.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.ticks.y.right = element_blank()
      )
    ),
    ggsignif.args = list(
      textsize = 2.5, tip_length = 0.01, color = "black",
      vjust = 0.5, size = 0.2
    )
  )

  # Save each plot to file
  ggsave(
    filename = paste0("metrics_Spectronaut/", score, "_FDR", FDR, ".png"),
    plot = p,
    dpi = 300,
    width = 6,
    height = 6
  )
}
