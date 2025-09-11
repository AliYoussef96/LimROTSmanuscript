
# Proteomics Differential Expression Analysis Pipeline

This repository contains a comprehensive pipeline for analyzing proteomics data from either MaxQuant or FragPipe, implementing multiple statistical methods for differential expression analysis.

## Overview

The pipeline implements nine different statistical methods:
- LimROTS (Reproducibility-Optimized Test Statistic with limma-trend)
- ROTS (Reproducibility-Optimized Test Statistic)
- Limma (Linear Models for Microarray and RNA-Seq Data)
- SAM (Significance Analysis of Microarrays)
- t-test
- ANOVA
- MSstats
- DEqMS (Differential Expression analysis of quantitative Mass Spectrometry data)
- DEP (Differential Enrichment analysis of Proteomics data)

## Directory Structure

```
.
├── data/                      # Contains experimental design information
│   └── data.info.csv         # Experiment metadata and design information
├── Maxquant/                 # MaxQuant input data directory
│   ├── *_LFQ_Maxquant_dlfq_pro_intensity.tsv
│   └── *_Maxquant_design.tsv
├── FragPipe/                 # FragPipe input data directory
│   ├── *_LFQ_FragPipe_dlfq_pro_intensity.tsv
│   └── *_FragPipe_design.tsv
├── Maxquant_results/         # Results from MaxQuant analysis
├── FragPipe_results/         # Results from FragPipe analysis
├── metrics_Maxquant/         # Evaluation metrics for MaxQuant
├── metrics_FragPipe/         # Evaluation metrics for FragPipe
└── Plot_files/              # Generated plots
    ├── MaxQuant/
    └── FragPipe/
```

## Requirements

1. R version ≥ 4.0.0
2. Required R packages:

```r
# Install required packages
install.packages(c(
    "stringr",
    "imputeLCMD",
    "dplyr",
    "ggplot2",
    "ggstatsplot",
    "ggsignif"
))

# Install Bioconductor packages
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

# Install SAM from CRAN
install.packages("samr", repos = "http://cran.us.r-project.org")
```

## Input Data Requirements

### Data Organization

The dataset can be downloaded from: [Download Link To Be Added]

1. Input files should be organized in `Maxquant/` or `FragPipe/` directories at the same level as your analysis scripts. These directories can be created in your working directory, or you can specify their location in the pipeline configuration.

2. Required input files for each experiment:
   - Intensity file: `*_LFQ_[Maxquant/FragPipe]_dlfq_pro_intensity.tsv`
   - Design file: `*_[Maxquant/FragPipe]_design.tsv`

### File Format Requirements

1. Intensity file format:
   - Tab-separated values
   - First column: Protein identifiers
   - Subsequent columns: Sample intensities
   - Column names should match sample names in design file

2. Design file format:
   - Tab-separated values
   - Required columns:
     - sample_name: Sample identifiers
     - condition: Experimental conditions (single characters, e.g., 'A', 'B')
     - replicate: Replicate numbers

3. data.info.csv format:
   - Required columns:
     - Dataset: Experiment identifiers
     - ID: Unique experiment identifiers

## Running the Pipeline

1. First, ensure your data is organized as described above.

2. Set the analysis mode in `run_analysis_pipeline.r`:
   ```r
   ANALYSIS_MODE <- "MaxQuant"  # or "FragPipe"
   ```

3. Run the complete pipeline:
   ```r
   source("run_analysis_pipeline.r")
   ```

   Or run individual scripts in sequence:
   ```r
   # For MaxQuant analysis:
   source("Maxquant.run.r")
   source("Maxquant.DEP.r")
   source("Maxquant.DEqMS.r")
   source("Maxquant.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   
   # For FragPipe analysis:
   source("FragPipe.run.r")
   source("FragPipe.DEP.r")
   source("FragPipe.DEqMS.r")
   source("FragPipe.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   ```

## Output Files

1. Analysis Results:
   - Stored in `[Maxquant/FragPipe]_results/`
   - Format: RDS files containing analysis results for each method
   - Naming: `[Method]_[Experiment]_[Contrast].rds`

2. Evaluation Metrics:
   - Stored in `metrics_[Maxquant/FragPipe]/`
   - Contains performance metrics for each analysis method
   - Includes: F1-Score, nMCC, G-Mean, etc.

3. Plots:
   - Stored in `Plot_files/[MaxQuant/FragPipe]/`
   - Generated visualizations for each metric
   - Format: High-resolution JPG files

## Troubleshooting

Common issues and solutions:

1. Missing directories:
   - The pipeline automatically creates required directories
   - Ensure write permissions in the working directory

2. Missing data files:
   - Check that input files exist in the correct directories
   - Verify file naming follows the required format

3. Package installation issues:
   - Ensure R version ≥ 4.0.0
   - Install Bioconductor packages using BiocManager
   - Check system dependencies for all packages

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added]

## License

[License information to be added]

