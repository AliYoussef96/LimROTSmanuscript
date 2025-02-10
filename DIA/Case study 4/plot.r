###############################################
######################  
###############################################
library(stringr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)

metrics <- read.csv("metrics/FDR_0.05")
colnames(metrics) <- str_remove(colnames(metrics), ".rds")

metrics <- metrics[metrics$X %in% c("F1_Score" , "nMCC" , "G_Mean"  ,
                                    "Balanced_Accuracy" , "pauc"),]

#metrics <- metrics[-c(8,9),]


metrics.melt <- reshape2::melt(metrics)
metrics.melt$DEA <- metrics.melt$variable

library(ggsci)

ggplot(metrics.melt , aes(x = X, y = value, fill = DEA )) + 
     geom_bar(stat = "identity" , position=position_dodge()) + 
     theme_bw() + scale_fill_lancet() + xlab("")

metrics.melt <- metrics.melt[metrics.melt$DEA %in% c("LimROTS" , "Limma"),]

library(dplyr)
library(tidyr)


ggsave("metrics/FDR0.05.tiff", dpi = 300, width = 10, height = 5)

metrics.melt$value <- as.numeric(data_long$value)


data_wide <- metrics.melt %>%
    pivot_wider(names_from = X, values_from = value)

# Add max and min rows for radar chart scaling
data_radar <- rbind(
    Max = rep(1, ncol(data_wide) - 1), # Maximum values for scaling
    Min = rep(0, ncol(data_wide) - 1), # Minimum values for scaling
    as.matrix(data_wide[, -1]) # Exclude Tool column
)

rownames(data_radar) <- c("Max", "Min", data_wide$DEA)

data_radar <- as.data.frame(data_radar)

data_radar[] <- sapply(data_radar[,], as.numeric)

row.names(data_radar)[3] <- "LimROTS"
row.names(data_radar)[4] <- "Limma"

library(fmsb)

radarchart(as.data.frame(data_radar), 
           axistype = 1, 
           pcol = c("#2980b9", "#e74c3c"), # Colors for the tools
           pfcol = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), # Transparent fill
           plwd = 2, # Line width
           cglcol = "grey", # Grid line color
           cglty = 1, # Grid line type
           axislabcol = "black", # Axis label color
           caxislabels = seq(0, 1, 0.2), # Axis labels
           vlcex = 0.8) # Text size for labels


library(ggplot2)

metrics.melt$variable <- droplevels(metrics.melt$variable)

ggplot(metrics.melt, aes(x = X, y = value, group = DEA, color = DEA, fill = DEA)) +
    geom_polygon(alpha = 0.4, show.legend = TRUE) +
    geom_line(size = 1) +
    geom_point(size = 3) +
    coord_polar() +
    theme_minimal() +
    theme(axis.text.x = element_text(size = 12, color = "black"),
          legend.position = "top") +
    labs(title = "Radar Chart: LimROTS vs Limma", y = "Score", x = "")




