# DIA Proteomics Analysis Pipeline - Case Study 2

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
1. R version ≥ 4.5.0

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

Example datasets can be downloaded from: [https://zenodo.org/records/17102211]

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

