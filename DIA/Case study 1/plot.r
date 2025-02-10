###############################################
######################  
###############################################
library(stringr)
library(dplyr)
library(ggplot2)
library(ggstatsplot)

metrics <- read.csv("metrics_Spectronaut/FDR_0.01")
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
        axis.text.x = element_text(color = "black", face = "bold", size = 7),
        axis.text.y = element_text(color = "black", face = "bold",  size = 9),
        axis.text = element_text(color = "black", face = "bold",  size = 7),
        text = element_text(color = "black", face = "bold"), 
        axis.title.y.right = element_blank(), 
              axis.text.y.right = element_blank(), 
              axis.ticks.y.right = element_blank()
    )), 
    ggsignif.args = list(textsize = 2.5, tip_length = 0.01 , color = "black",
                         vjust = 0.5 , size = 0.2 )
) 

ggsave(paste0( "metrics_Spectronaut/" , score, "_FDR" , FDR , ".png"), plot = p, dpi = 300, width = 6, height = 6)

}
# 
# metrics.melt[is.na(metrics.melt)] <- 0
# 
# metrics.melt$exp <- as.factor(metrics.melt$exp)
# 
# p <- grouped_ggbetweenstats(metrics.melt,
#                        x = DEA,
#                        y = value,
#                        grouping.var = exp,
#                        type = "nonparametric",
#                        centrality.type = "nonparametric" ,
#                        p.adjust.method = "BH",
#                        pairwise.display = "s",
#                        package = "ggsci",
#                        palette = "lanonc_lancet")
# 
# ggsave(paste0( "metrics_Spectronaut/", "ALL_FDR" , FDR , ".tiff"), plot = p, dpi = 300, width = 20, height = 20)

############# ###################### ##########
############## PCA

#metrics[is.na(metrics)] <- 0

row.names(metrics) <- metrics$X

metrics <- metrics[,-1]
metrics <- t(metrics)
metrics <- data.frame(metrics)
metrics$DEA <- str_split_fixed(row.names(metrics) , fixed("_"),2)[,1]
metrics <- na.omit(metrics)

pca_result <- prcomp(metrics[-6], scale. = TRUE)
plot(pca_result)
summary(pca_result)
pca_scores <- as.data.frame(pca_result$x)
pca_scores$DEA <- metrics$DEA

ggplot(pca_scores, aes(x = PC1, y = PC2, color = DEA)) +
    geom_point(size = 3) +
    ggtitle("PCA of Statistical Methods") +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    theme_minimal()

loadings <- as.data.frame(pca_result$rotation)


# Perform k-means clustering (2 clusters for good/bad methods)
set.seed(42)
clusters <- kmeans(pca_scores[, c("PC1", "PC2")], centers = 3)

# Add cluster assignments to the PCA scores
pca_scores$Cluster <- as.factor(clusters$cluster)
pca_scores$DEA <- metrics$DEA

# Plot PCA with clusters
ggplot(pca_scores, aes(x = PC1, y = PC2)) +
    geom_point(size = 3, aes(color = DEA, shape = Cluster)) +
    stat_ellipse(aes(fill = Cluster), type = "norm", alpha = 0.2, geom = "polygon") +
    ggtitle("PCA of Statistical Methods with Clusters") +
    xlab("Principal Component 1") +
    ylab("Principal Component 2") +
    theme_minimal()

############# ###################### ##########

############ MDS


fold_change_distance <- as.matrix(dist(log2(metrics[-6] + 1)))  # Adding 1 to avoid log(0)
# Perform MDS
mds_result <- cmdscale(fold_change_distance, k = 2)

# Convert MDS results to a dataframe
mds_data <- as.data.frame(mds_result)
colnames(mds_data) <- c("Dim1", "Dim2")
mds_data$Method <- rownames(metrics)

# Perform k-means clustering (2 clusters)
set.seed(42)
clusters <- kmeans(mds_data[, c("Dim1", "Dim2")], centers = 3)

# Add cluster assignments
mds_data$Cluster <- as.factor(clusters$cluster)
mds_data$DEA <- metrics$DEA

# MDS plot with ellipses around clusters
ggplot(mds_data, aes(x = Dim1, y = Dim2)) +
    geom_point(size = 3, aes(colour = DEA, shape = Cluster) ) +
    stat_ellipse(aes(fill = Cluster), type = "t", alpha = 0.2, geom = "polygon") +
    ggtitle("MDS of Methods Based on Fold-Change Distance") +
    xlab("Dimension 1") +
    ylab("Dimension 2") +
    theme_minimal()

############# ###################### ##########
############# ###################### ##########
############# ###################### ##########

metrics.melt$comp <- str_split_fixed(metrics.melt$variable , fixed("_"),3)[,3]
metrics.melt$comp <- ifelse( grepl("_", metrics.melt$comp , fixed = TRUE) , str_split_fixed(metrics.melt$comp  , fixed("_"),2)[,2], metrics.melt$comp)

which.max.dea <- data.frame()
which.min.dea <- data.frame()

for (i in unique(metrics.melt$exp)){
    
    d.temp <- metrics.melt[metrics.melt$exp == i, ]
    
    for(q in unique(d.temp$comp)){
        d.temp.temp <- d.temp[d.temp$comp == q,]
        
    for(j in unique(d.temp.temp$X)){
        d.temp.temp.temp <- d.temp.temp[d.temp.temp$X == j,]
        
        max.dea <- as.vector(d.temp.temp.temp[which.max(d.temp.temp.temp$value),]$DEA)
        min.dia <- as.character(d.temp.temp.temp[which.min(d.temp.temp.temp$value),]$DEA)
        
        which.max.dea <- rbind(which.max.dea, data.frame(exp = i,
                                                         comp = q,
                                                         index = j,
                                                         
                                                         max.dea = max.dea))
        
        which.min.dea <- rbind(which.min.dea, data.frame(exp = i,
                                                         comp = q,
                                                         index = j,
                                                         min.dia = min.dia))
    }
        
    }
}

table(which.max.dea$max.dea)
table(which.min.dea$min.dia)


table(which.max.dea$index , which.max.dea$max.dea)

table(which.min.dea$index , which.min.dea$min.dia)
