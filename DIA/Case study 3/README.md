# DIA Proteomics Analysis Pipeline - Case Study 3


## Overview

The pipeline implements multiple statistical methods for differential expression analysis:
- LimROTS (Reproducibility-Optimized Test Statistic with limma-trend)
- ROTS (Reproducibility-Optimized Test Statistic)
- Limma (Linear Models for Microarray and RNA-Seq Data)
- MSstats (with DIA-specific processing)
- DEqMS (Differential Expression analysis of quantitative Mass Spectrometry data)

## Directory Structure

```
.
├── data/                      # Contains experimental design information
│   └── data.info.csv         # Experiment metadata and design information
├── DIANN/                    # DIA-NN input data directory
│   ├── HEof_DIANN_dlfq.tsv       # Quantification data
│   ├── HEof_DIANN_design.tsv     # Design file
│   └── HEof_DIANN_design_msstats.tsv  # MSstats annotations
├── Spectronaut/             # Spectronaut input data directory
│   ├── HEof_spt_dlfq.tsv        # Quantification data
│   ├── HEof_spt_design.tsv      # Design file
│   └── HEof_spt_design_msstats.tsv   # MSstats annotations
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
    "reshape2"
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
    "limma"
))
```

## Input Data Requirements

### Data Organization

Place your input files in the appropriate directories:

1. DIA-NN data:
   - Place files in `/DIANN/` directory:
     * `HEof_DIANN_dlfq.tsv`: Quantification data
     * `HEof_DIANN_design.tsv`: Experimental design
     * `HEof_DIANN_design_msstats.tsv`: MSstats annotations

2. Spectronaut data:
   - Place files in `/Spectronaut/` directory:
     * `HEof_spt_dlfq.tsv`: Quantification data
     * `HEof_spt_design.tsv`: Experimental design
     * `HEof_spt_design_msstats.tsv`: MSstats annotations

### File Format Requirements

1. Quantification data files (dlfq.tsv):
   - Tab-separated values
   - Protein identifiers in first column
   - Sample intensities in subsequent columns

2. Design files:
   - Tab-separated values
   - Must include columns: Run, Condition, BioReplicate

3. MSstats annotation files:
   - Tab-separated values
   - Must follow MSstats format requirements

## Running the Pipeline

1. Ensure all required packages are installed
2. Place input files in appropriate directories
3. Execute the analysis pipeline:
```r
source("run_analysis_pipeline.r")
```

## Output

The pipeline generates:
1. Statistical analysis results in results directories
2. Performance metrics in metrics directories
3. Visualization plots in Plot_files directory

## Performance Metrics

The pipeline evaluates:
- F1 Score
- Matthews Correlation Coefficient (MCC)
- G-Mean
- Balanced Accuracy
- Partial Area Under Curve (pAUC)
- Fowlkes-Mallows Index

