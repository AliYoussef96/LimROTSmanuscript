# LimROTS Manuscript: Reproducible Analysis Guide

## Quick Start Guide

1. Clone this repository
2. Install required R packages (see Dependencies section below)
3. Run your chosen analysis:
   - DDA Analysis
   - DIA Analysis (Case Studies 1-4)
   - UPenn Cohort Analysis

Note: All the required data files and code, included in their respective directories, are available for direct download from https://zenodo.org/records/17102211.

## Dependencies

```R
# Install required R packages
install.packages(c("stringr", "ggplot2", "dplyr" , "samr"))
BiocManager::install(c("DEP", "limma" ,"DEqMS", "MSstats", "LimROTS" ,"SummarizedExperiment"))
```

## Overview

This repository provides all necessary R scripts to reproduce the analyses presented in the LimROTS manuscript.

### Methods Compared

- LimROTS
- ROTS
- Limma
- SAM
- t-test
- ANOVA
- MSstats
- DEqMS
- DEP

---

## DDA Analysis

The `DDA/` directory contains R scripts for reproducing all results from the manuscript's DDA-based proteomics analysis. See `DDA/readme.md` for detailed information about the DDA analysis workflow.

### Running the Analysis

The simplest way to run the complete DDA analysis is to use the automated pipeline:
```R
cd DDA
Rscript run_analysis_pipeline.r
```

This script will automatically:
1. Run all benchmark methods (FragPipe and MaxQuant)
2. Execute all statistical analyses
3. Perform evaluations
4. Generate result metrics

### Output Structure

- Results will be stored in:
  - `FragPipe_results/` and `Maxquant_results/` for analysis outputs
  - `metrics_FragPipe/` and `metrics_Maxquant/` for evaluation metrics
  (All folders are created automatically)

---

## DIA Analysis

The `DIA/` directory includes scripts to reproduce the DIA-based analyses from Case Studies 1 through 4. Each case study has its own detailed documentation:
- Case Study 1: See `DIA/Case study 1/README.MD`
- Case Study 2: See `DIA/Case study 2/README.MD`
- Case Study 3: See `DIA/Case study 3/README.md`
- Case Study 4: See `DIA/Case study 4/README.md`

### Running the Analyses

Each case study can be run using its automated pipeline script:

```R
# Case Study 1
cd "DIA/Case study 1"
Rscript run_analysis_pipeline.r

# Case Study 2
cd "DIA/Case study 2"
Rscript run_analysis_pipeline.r

# Case Study 3
cd "DIA/Case study 3"
Rscript run_analysis_pipeline.r

# Case Study 4
cd "DIA/Case study 4"
Rscript run_analysis_pipeline.r
```

Each `run_analysis_pipeline.r` script will automatically:
1. Process both DIA-NN and Spectronaut data
2. Run all statistical methods (MSstats, DEqMS, DEP)
3. Perform evaluations
4. Generate result metrics

### Output Structure

Results for each case study will be stored in:
- `spectronaut_results/` and `DIANN_results/` for analysis outputs
- `metrics_spectronaut/` and `metrics_DIANN/` for evaluation metrics
(All folders are created automatically)

---

## UPenn Cohort Analysis

### Dataset Access

- **Source:** [Synapse: syn20933797](https://www.synapse.org/Synapse:syn20933797/wiki/596247)

### Script Directory

- **Folder:** `ConsensusProteinCoexpressionStudyPenn/`
- **Purpose:** Contains preprocessing and analysis code for the UPenn Alzheimer’s disease cohort.
- **Methods Used:** LimROTS, Limma, and ROTS.

---

## Notes

- All output folders are created automatically by the scripts
- All R scripts assume that the required dependencies and packages are installed

---

## Citation

If you use this repository, please cite the LimROTS manuscript accordingly.

