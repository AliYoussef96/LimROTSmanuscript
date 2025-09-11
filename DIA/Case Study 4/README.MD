# DIA Proteomics Analysis Pipeline - Case Study 4

## Overview

The pipeline implements statistical methods optimized for HEof data analysis:
- LimROTS (Reproducibility-Optimized Test Statistic with limma-trend)
- Limma (Linear Models for Microarray and RNA-Seq Data)
- DEqMS (Differential Expression analysis of quantitative Mass Spectrometry data)

## Directory Structure

```
.
├── results/                  # Analysis results directory
├── metrics/                  # Performance metrics
└── Plot_files/              # Generated plots and visualizations
```

Note: This case study uses the built-in UPS1.Case4 dataset which is loaded directly using `data(UPS1.Case4)`. No input data files are required.

## Requirements

### Software Requirements
1. R version ≥ 4.0.0
2. Required R packages (see below)

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
    "BiocParallel",
    "LimROTS",
    "DEqMS",
    "limma"
))
```

## Input Data Requirements

### Data Organization

This case study uses the built-in UPS1.Case4 dataset which contains:

1. Expression Data:
   - Protein quantification matrix
   - Sample measurements across conditions
   - Pre-processed and normalized data

2. Sample Metadata:
   - Concentration information (Conc.)
   - Tool information
   - Batch information (fake.batch)

3. Feature Metadata:
   - Protein/Gene identifiers
   - Annotation information

The data is accessed directly through R using:
```r
data(UPS1.Case4)
```

### Data Structure

The UPS1.Case4 dataset is a SummarizedExperiment object containing:
1. Assay data (expression matrix)
2. Column data (sample information)
3. Row data (feature annotations)

## Running the Pipeline

1. Install required packages:
```r
source("install_dependencies.r")  # If provided
```

2. Place input files in the data directory

3. Execute the analysis pipeline:
```r
source("run_analysis_pipeline.r")
```

## Pipeline Steps

1. HEof Analysis (HEof.r):
   - Loads and processes input data
   - Performs differential expression analysis
   - Generates initial results

2. Evaluation (evaluate.r):
   - Calculates performance metrics
   - Assesses statistical significance
   - Generates evaluation reports

3. Visualization (plot.r):
   - Creates performance plots
   - Generates comparison visualizations
   - Outputs summary figures

## Output

The pipeline generates three types of output:

1. Results (in /results):
   - Differential expression analysis results
   - Statistical test outcomes
   - Significance measurements

2. Metrics (in /metrics):
   - F1 Score
   - Matthews Correlation Coefficient (MCC)
   - G-Mean
   - Balanced Accuracy
   - Partial Area Under Curve (pAUC)
   - Fowlkes-Mallows Index

3. Plots (in /Plot_files):
   - Performance metric visualizations
   - Comparison plots
   - Summary figures

## Performance Metrics

The pipeline evaluates method performance using:
- F1 Score: Balance between precision and recall
- nMCC: Normalized Matthews Correlation Coefficient
- G-Mean: Geometric mean of sensitivity and specificity
- Balanced Accuracy: Average of sensitivity and specificity
- pAUC: Partial Area Under the ROC Curve
- Fowlkes-Mallows Index: Geometric mean of precision and recall

