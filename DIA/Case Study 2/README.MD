# DIA Proteomics Analysis Pipeline - Case Study 2

This repository contains a comprehensive pipeline for analyzing DIA proteomics data from either DIA-NN or Spectronaut, implementing multiple statistical methods for differential expression analysis with a focus on batch effect handling.

## Overview

The pipeline implements multiple statistical methods optimized for DIA data with batch effect consideration:
- LimROTS (Reproducibility-Optimized Test Statistic with limma-trend)
- ROTS (Reproducibility-Optimized Test Statistic)
- Limma (Linear Models for Microarray and RNA-Seq Data)
- SAM (Significance Analysis of Microarrays)
- t-test
- ANOVA
- MSstats (with DIA-specific processing)
- DEqMS (Differential Expression analysis of quantitative Mass Spectrometry data)
- DEP (Differential Enrichment analysis of Proteomics data)

## Directory Structure

```
.
├── data/                      # Contains experimental design information
│   └── data.info.csv         # Experiment metadata and design information
├── DIANN/                    # DIA-NN input data directory
│   ├── HEof_n600_DIA_DIANN_dlfq.tsv    # Narrow window data
│   ├── HEof_w600_DIA_DIANN_dlfq.tsv    # Wide window data
│   ├── HEof_*_DIA_DIANN_design.tsv     # Design files
│   └── HEof_*_DIA_DIANN_design_msstats.tsv  # MSstats annotations
├── Spectronaut/             # Spectronaut input data directory
│   ├── HEof_n600_DIA_spt_dlfq.tsv     # Narrow window data
│   ├── HEof_w600_DIA_spt_dlfq.tsv     # Wide window data
│   ├── HEof_*_DIA_spt_design.tsv      # Design files
│   └── HEof_*_DIA_spt_design_msstats.tsv   # MSstats annotations
├── DIANN_results/          # Results from DIA-NN analysis
├── Spectronaut_results/    # Results from Spectronaut analysis
├── metrics_DIANN/          # Evaluation metrics for DIA-NN
├── metrics_Spectronaut/    # Evaluation metrics for Spectronaut
└── Plot_files/             # Generated plots
    ├── DIANN/              # DIA-NN specific plots
    └── Spectronaut/        # Spectronaut specific plots
```

## Requirements

### Software Requirements
1. R version ≥ 4.0.0
2. DIA-NN (≥ 1.8) or Spectronaut (≥ 15.0) for raw data processing

### R Package Dependencies

```r

# Core dependencies
install.packages(c(
    "stringr",
    "data.table",
    "tidyr",
    "dplyr",
    "ggplot2",
    "ggstatsplot",
    "ggsignif",
    "ggrepel",
    "pheatmap",
    "RColorBrewer"
))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
    "SummarizedExperiment",
    "ROTS",
    "BiocParallel",
    "LimROTS",
    "MSstats",
    "DEqMS",
    "DEP",
    "limma"
))

# Additional dependencies
install.packages(c(
    "samr",
    "imputeLCMD",
    "rrcovNA"
), repos = "http://cran.us.r-project.org")
```

## Input Data Requirements

### Data Organization

Example datasets can be downloaded from: [Download Link To Be Added]

1. DIA-NN data requirements:
   - Place DIA-NN quantification files in `/DIANN/` directory:
     * `HEof_n600_DIA_DIANN_dlfq.tsv`: Narrow window data
     * `HEof_w600_DIA_DIANN_dlfq.tsv`: Wide window data
   - Design files in `/DIANN/` directory:
     * `HEof_*_DIA_DIANN_design.tsv`: Experimental design
     * `HEof_*_DIA_DIANN_design_msstats.tsv`: MSstats annotations

2. Spectronaut data requirements:
   - Place Spectronaut quantification files in `/Spectronaut/` directory:
     * `HEof_n600_DIA_spt_dlfq.tsv`: Narrow window data
     * `HEof_w600_DIA_spt_dlfq.tsv`: Wide window data
   - Design files in `/Spectronaut/` directory:
     * `HEof_*_DIA_spt_design.tsv`: Experimental design
     * `HEof_*_DIA_spt_design_msstats.tsv`: MSstats annotations

## Analysis Scripts

The analysis pipeline consists of several R scripts:

1. `DIANN.DEqMS.r` - Performs DEqMS analysis on DIA-NN data
2. `DIANN.MSstats.r` - Runs MSstats analysis on DIA-NN data
3. `Spectronaut.DEqMS.r` - Performs DEqMS analysis on Spectronaut data
4. `Spectronaut.MSstats.r` - Runs MSstats analysis on Spectronaut data
5. `evaluate.r` - Evaluates performance metrics
6. `plot.r` - Generates performance comparison plots

## Results

The analysis generates several types of outputs:

1. Statistical analysis results (in `DIANN_results/` or `Spectronaut_results/`):
   - DEqMS analysis results
   - MSstats analysis results

2. Performance metrics (in `metrics_DIANN/` or `metrics_Spectronaut/`):
   - F1 Score
   - nMCC (normalized Matthews Correlation Coefficient)
   - G-Mean
   - Balanced Accuracy
   - pAUC (partial Area Under Curve)

3. Visualization (in `Plot_files/`):
   - Performance comparison plots
   - Method comparison plots
   - Batch effect visualization

## Example Plots

Below are example plots generated from the analysis showing the performance metrics for different methods:


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

![](Plot_files/Balanced_Accuracy_FDR0.05.jpg)


![](Plot_files/F1_Score_FDR0.05.jpg)


![](Plot_files/G_Mean_FDR0.05.jpg)


![](Plot_files/nMCC_FDR0.05.jpg)


![](Plot_files/pauc_FDR0.05.jpg)
