# DIA Proteomics Analysis Pipeline DIA Case Study 1

## Overview

The pipeline implements nine different statistical methods specifically optimized for DIA data:
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
│   ├── report.tsv            # DIA-NN output files
│   ├── *_DIANN_design.tsv   # Experimental design files
│   └── *_DIANN_design_msstats.tsv  # MSstats annotation files
├── Spectronaut/             # Spectronaut input data directory
│   ├── Report.tsv           # Spectronaut output files
│   ├── *_spt_design.tsv    # Experimental design files
│   └── *_spt_design_msstats.tsv   # MSstats annotation files
├── DIANN_results/          # Results from DIA-NN analysis
├── Spectronaut_results/    # Results from Spectronaut analysis
├── metrics_DIANN/          # Evaluation metrics for DIA-NN
├── metrics_Spectronaut/    # Evaluation metrics for Spectronaut
└── Plot_files/             # Generated plots
    ├── DIANN/
    └── Spectronaut/
```


## Requirements

### Software Requirements
1. R version ≥ 4.5

### R Package Dependencies

```r
# Core dependencies
install.packages(c(
    "stringr",
    "data.table",
    "tidyr",
    "dplyr",
    "ggplot2",
    "patchwork",
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

Example datasets can be downloaded from: https://zenodo.org/records/17102211

1. DIA-NN data requirements:
   - Place DIA-NN report.tsv files in `/DIANN/{experiment}/` directories
   - Design files in `/DIANN/` directory:
     * `{experiment}_DIANN_design.tsv`: Experimental design
     * `{experiment}_DIANN_design_msstats.tsv`: MSstats annotations

2. Spectronaut data requirements:
   - Place Spectronaut Report.tsv files in `/Spectronaut/{experiment}/` directories
   - Design files in `/Spectronaut/` directory:
     * `{experiment}_spt_design.tsv`: Experimental design
     * `{experiment}_spt_design_msstats.tsv`: MSstats annotations

### File Format Requirements

1. DIA-NN/Spectronaut report format:
   - Tab-separated values
   - Required columns:
     * Protein identifiers (PG.ProteinAccessions)
     * Precursor identifiers (EG.PrecursorId)
     * Quantitative values (PG.Quantity)
     * Run information (R.FileName, R.Condition, R.Replicate)

2. Design file format:
   - Tab-separated values
   - Required columns:
     * sample_name: Sample identifiers
     * condition: Experimental conditions (single characters, e.g., 'A', 'B')
     * replicate: Replicate numbers

3. MSstats design format:
   - Tab-separated values
   - Additional MSstats-specific columns
   - Matches MSstats input requirements for DIA data

## Running the Pipeline

1. First, ensure your data is organized as described above.

2. Set the analysis mode in `run_analysis_pipeline.r`:
   ```r
   ANALYSIS_MODE <- "DIANN"  # or "Spectronaut"
   ```

3. Run the complete pipeline:
   ```r
   source("run_analysis_pipeline.r")
   ```

   Or run individual scripts in sequence:
   ```r
   # For DIA-NN analysis:
   source("DIANN.run.r")
   source("DIANN.DEP.r")
   source("DIANN.DEqMS.r")
   source("DIANN.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   
   # For Spectronaut analysis:
   source("spectronaut.run.r")
   source("spectronaut.DEP.r")
   source("spectronaut.DEqMS.r")
   source("spectronaut.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   ```

## Output Files

1. Analysis Results:
   - Stored in `[DIANN/Spectronaut]_results/`
   - Format: RDS files containing analysis results for each method
   - Naming: `[Method]_[Experiment]_[Contrast].rds`

2. Raw MSstats Files:
   - Stored in `rawMSstats/`
   - Contains formatted input data for MSstats analysis
   - Includes protein and peptide-level quantification

3. Evaluation Metrics:
   - Stored in `metrics_[DIANN/Spectronaut]/`
   - Contains performance metrics for each analysis method
   - Includes: F1-Score, nMCC, G-Mean, etc.

4. Plots:
   - Stored in `Plot_files/[DIANN/Spectronaut]/`
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
   - Ensure DIA-NN/Spectronaut exports are in the correct format

3. Package installation issues:
   - Ensure R version ≥ 4.0.0
   - Install Bioconductor packages using BiocManager
   - Check system dependencies for all packages

4. MSstats-specific issues:
   - Verify MSstats annotation files match sample names exactly
   - Check for missing values in required columns
   - Ensure proper formatting of condition and replicate information

## Citation

If you use this pipeline in your research, please cite:
[Citation information to be added]

## License

[License information to be added]
   - Required columns:
     - Dataset: Experiment identifiers
     - ID: Unique experiment identifiers

## Running the Pipeline

1. First, ensure your data is organized as described above.

2. Set the analysis mode in `run_analysis_pipeline.r`:
   ```r
   ANALYSIS_MODE <- "DIANN"  # or "Spectronaut"
   ```

3. Run the complete pipeline:
   ```r
   source("run_analysis_pipeline.r")
   ```

   Or run individual scripts in sequence:
   ```r
   # For DIA-NN analysis:
   source("DIANN.run.r")
   source("DIANN.DEP.r")
   source("DIANN.DEqMS.r")
   source("DIANN.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   
   # For Spectronaut analysis:
   source("spectronaut.run.r")
   source("Spectronaut.DEP.r")
   source("Spectronaut.DEqMS.r")
   source("Spectronaut.MSstats.r")
   source("evaluate.r")
   source("plot.r")
   ```

## Output Files

1. Analysis Results:
   - Stored in `[DIANN/Spectronaut]_results/`
   - Format: RDS files containing analysis results for each method
   - Naming: `[Method]_[Experiment]_[Contrast].rds`

2. Evaluation Metrics:
   - Stored in `metrics_[DIANN/Spectronaut]/`
   - Contains performance metrics for each analysis method
   - Includes: F1-Score, nMCC, G-Mean, etc.

3. Plots:
   - Stored in `Plot_files/[DIANN/Spectronaut]/`
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


