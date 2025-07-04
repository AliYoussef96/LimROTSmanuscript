```r
    library(stringr)
    library(dplyr)
    library(ggplot2)
    library(ggstatsplot)
    library(ggsignif)

    metrics <- read.csv("metrics_Spectronaut//FDR_0.05.csv")
    colnames(metrics) <- str_remove(colnames(metrics), ".rds")
    colnames(metrics) <- str_remove(colnames(metrics), "_LFQ")
    metrics <- metrics[metrics$X %in% c("F1_Score" , "nMCC" , "G_Mean" , "Fowlkes_Mallows_Index" ,
                                        "Balanced_Accuracy" , "pauc"),]


    metrics.melt <- reshape2::melt(metrics)

    ## Using X as id variables

    metrics.melt$DEA <- str_split_fixed(metrics.melt$variable , fixed("_"),2)[,1]
    metrics.melt$exp <- str_split_fixed(metrics.melt$variable , fixed("_"),2)[,2]
    metrics.melt$exp <- sub("_.*", "", metrics.melt$exp)

    metrics.melt <- metrics.melt[metrics.melt$DEA %in% c("DEP", "Limma", "LimROTS", "ANOVA"),]

    FDR <- 0.05

    for(score in unique(metrics.melt$X)){
        
        df <-  metrics.melt[metrics.melt$X == score,]
        
        df[is.na(df)] <- 0
        
        p <- ggbetweenstats(
            data = df,
            x = DEA,
            y = value,
            type = "robust",
            centrality.type = "nonparametric" ,
            p.adjust.method = "BH",
            pairwise.display = "none",
            bf.message = FALSE,
            results.subtitle = FALSE,
            package = "ggsci",
            palette = "lanonc_lancet",
            ylab = score,
            centrality.label.args = list(size = 4.5, segment.linetype = 1,
                                         nudge_x = 0.25),
            ggplot.component = list(  theme(
                axis.text.x = element_text(color = "black", face = "bold", size = 12),
                axis.text.y = element_text(color = "black", face = "bold",  size = 12),
                axis.text = element_text(color = "black", face = "bold",  size = 12),
                axis.title.y.left  = element_text(color = "black", face = "bold", size = 17),
                axis.title.x.bottom   = element_text(color = "black", face = "bold", size = 17),
                text = element_text(color = "black", face = "bold"),
                legend.title = element_text(color = "black", face = "bold", size = 13.5),
                legend.text = element_text(color = "black", face = "bold", size = 13.5),
                axis.title.y.right = element_blank(), 
                axis.text.y.right = element_blank(), 
                axis.ticks.y.right = element_blank()
            )), 
            ggsignif.args = list(textsize = 25, tip_length = 0.01 , color = "black",
                                 vjust = 0.5 , size = 0.2 )
        )  + scale_x_discrete(labels = levels(as.factor(df$DEA))) +
          ggplot2::scale_color_manual(values = c("#00468BFF", "#ED0000FF" , "#0099B4FF", "#925E9FFF"))
        
        print(p)

        
    }
```

![](Plot_files/figure-markdown_strict/unnamed-chunk-1-1.png)


![](Plot_files/figure-markdown_strict/unnamed-chunk-1-2.png)


![](Plot_files/figure-markdown_strict/unnamed-chunk-1-3.png)


![](Plot_files/figure-markdown_strict/unnamed-chunk-1-4.png)


![](Plot_files/figure-markdown_strict/unnamed-chunk-1-5.png)


![](Plot_files/figure-markdown_strict/unnamed-chunk-1-6.png)
