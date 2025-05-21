# Load required libraries for string manipulation, data manipulation, and plotting
library(stringr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)
library(ggradar)

# Read the CSV file containing metric data with FDR 0.05 threshold
metrics <- read.csv("metrics/FDR_0.05")

# Remove the substring ".rds" from column names in the metrics dataframe
colnames(metrics) <- str_remove(colnames(metrics), ".rds")

# Subset the metrics dataframe to include only specific performance metrics of interest
metrics <- metrics[metrics$X %in% c("F1_Score" , "nMCC" , "G_Mean"  ,
                                    "Balanced_Accuracy" , "pauc"),]


# Reshape the metrics dataframe from wide to long format for plotting
metrics.melt <- reshape2::melt(metrics)

# Create a new column 'DEA' to store the method/tool names from the variable column
metrics.melt$DEA <- metrics.melt$variable

# Load color palette library for ggplot
library(ggsci)

# Create a bar plot comparing performance metrics across different methods/tools
ggplot(metrics.melt , aes(x = X, y = value, fill = DEA )) + 
     geom_bar(stat = "identity" , position=position_dodge()) + 
     theme_bw() + scale_fill_lancet() + xlab("")

# Filter the melted data to include only two specific methods: LimROTS and Limma
metrics.melt <- metrics.melt[metrics.melt$DEA %in% c("LimROTS" , "Limma"),]

# Load libraries for further data manipulation
library(dplyr)
library(tidyr)

# Save the current plot as a high-resolution TIFF image
ggsave("metrics/FDR0.05.tiff", dpi = 300, width = 10, height = 5)

# Convert the 'value' column in metrics.melt to numeric to ensure proper numeric operations
metrics.melt$value <- as.numeric(data_long$value)

# Pivot the data from long format to wide format with metrics as columns
data_wide <- metrics.melt %>%
    pivot_wider(names_from = X, values_from = value)

# Prepare data for radar chart plotting by adding max and min rows for scaling purposes
data_radar <- rbind(
    Max = rep(1, ncol(data_wide) - 1), # Maximum values set to 1 for each metric
    Min = rep(0, ncol(data_wide) - 1), # Minimum values set to 0 for each metric
    as.matrix(data_wide[, -1]) # Actual metric values, excluding the first column (DEA/tool names)
)

# Assign row names for radar chart data: Max, Min, and tool names
rownames(data_radar) <- c("Max", "Min", data_wide$DEA)

# Convert the data frame to numeric to avoid issues during plotting
data_radar <- as.data.frame(data_radar)
data_radar[] <- sapply(data_radar[,], as.numeric)

# Rename row names explicitly to match the tools LimROTS and Limma
row.names(data_radar)[3] <- "LimROTS"
row.names(data_radar)[4] <- "Limma"

# Load the fmsb package used for radar chart plotting
library(fmsb)

# Generate a radar chart with specified aesthetics and customization
radarchart(as.data.frame(data_radar), 
           axistype = 1, 
           pcol = c("#2980b9", "#e74c3c"), # Colors for the lines representing each tool
           pfcol = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), # Transparent fill colors for polygons
           plwd = 2, # Line width of polygons
           cglcol = "grey", # Color of the circular grid lines
           cglty = 1, # Type of the circular grid lines (solid)
           axislabcol = "black", # Color of the axis labels
           caxislabels = seq(0, 1, 0.2), # Labels on the axis ranging from 0 to 1
           vlcex = 0.8) # Size of the variable labels around the radar plot


# Drop unused factor levels in the 'variable' column to clean data
metrics.melt$variable <- droplevels(metrics.melt$variable)

# Create a radar chart using ggplot2 with polar coordinates, comparing LimROTS and Limma
ggplot(metrics.melt, aes(x = X, y = value, group = DEA, fill = DEA)) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    coord_polar() + # Convert plot to polar coordinates (radar chart)
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          legend.position = "top") +
    labs(title = "Radar Chart: LimROTS vs Limma", y = "Score", x = "")

# Remove Max and Min rows from radar chart data for ggradar compatibility
data_radar <- data_radar[-c(1,2),]

# Assign tool names as a new first column for ggradar input
data_radar[,1] <- c("LimROTS","Limma")

# Reorder rows to match desired order (Limma first, LimROTS second)
data_radar <- data_radar[c(2, 1), ]

# Display the final data_radar table in console
data_radar

# Generate radar chart using ggradar with customized point sizes and grid labels
p <- ggradar(data_radar,    
             grid.label.size = 7,  # Size of grid annotation labels (e.g., 0%, 50%)
             axis.label.size = 5,  # Size of variable names on axes
             group.point.size = 5, # Size of points representing groups
             values.radar = seq(0,1,0.5) )  # Values shown on radar grid

# Display the radar chart plot
p

# Save the radar chart plot as a TIFF file with specified resolution and size
ggsave("metrics//radar.plot.tiff",p , dpi = 300, height = 10, width = 10)
