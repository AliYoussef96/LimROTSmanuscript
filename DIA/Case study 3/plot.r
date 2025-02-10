###############################################
######################  
###############################################
library(stringr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)

metrics <- read.csv("metrics_DIANN/FDR_0.01")
colnames(metrics) <- str_remove(colnames(metrics), ".rds")

metrics <- metrics[metrics$X %in% c("F1_Score" , "nMCC" , "G_Mean" , "Fowlkes_Mallows_Index" ,
                                    "Balanced_Accuracy" , "pauc"),]

#metrics <- metrics[-c(8,9),]


metrics.melt <- reshape2::melt(metrics)
metrics.melt$DEA <- str_split_fixed(metrics.melt$variable , fixed("_"),2)[,1]

metrics.melt$exp <- str_split_fixed(metrics.melt$variable , fixed("_"),2)[,2]
metrics.melt$exp <- sub("_.*", "", metrics.melt$exp)


median.DEA <- metrics.melt %>% group_by(X, DEA) %>% summarise(median.DEA = median(value , na.rm = TRUE))
median.DEA <- data.frame(median.DEA)

# ggplot(median.DEA , aes(x = median.DEA, y = X, fill = DEA )) + 
#     geom_bar(stat = "identity" , position=position_dodge()) + 
#     theme_bw()



median.DEA <- metrics.melt %>% group_by(X, DEA, exp) %>% summarise(median.DEA = median(value , na.rm = TRUE))
median.DEA <- data.frame(median.DEA)

# ggplot(median.DEA , aes(x = median.DEA, y = X, fill = DEA)) + 
#     geom_bar(stat = "identity" , position=position_dodge()) + 
#     theme_bw() + facet_grid(.~exp)




FDR <- 0.01

#####################
#####################



for(score in unique(metrics.melt$X)){

df <-  metrics.melt[metrics.melt$X == score,]

df[is.na(df)] <- 0

p <- ggbetweenstats(
    data = df,
    x = DEA,
    y = value,
    type = "nonparametric",
    centrality.type = "nonparametric" ,
    p.adjust.method = "BH",
    pairwise.display = "s",
    bf.message = FALSE,
    results.subtitle = FALSE,
    package = "ggsci",
    palette = "lanonc_lancet",
    ylab = score,
    ggplot.component = list(  theme(
        axis.text.x = element_text(color = "black", face = "bold", size = 10),
        axis.text.y = element_text(color = "black", face = "bold",  size = 10),
        axis.text = element_text(color = "black", face = "bold",  size = 10),
        text = element_text(color = "black", face = "bold"), 
        axis.title.y.right = element_blank(), 
        axis.text.y.right = element_blank(), 
        axis.ticks.y.right = element_blank()
    )), 
    ggsignif.args = list(textsize = 3, tip_length = 0.01 , color = "black",
                         vjust = 0.5 , size = 0.2 )
) 

ggsave(paste0( "metrics_DIANN/" , score, "_FDR" , FDR , ".png"), plot = p, dpi = 300, width = 6, height = 6)
}

