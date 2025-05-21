# LimROTS Manuscript: Reproducible Analysis Guide

## Overview

This repository provides all necessary R scripts to reproduce the analyses presented in the LimROTS manuscript. It includes a comprehensive comparison of differential expression methods applied to both Data-Dependent Acquisition (DDA) and Data-Independent Acquisition (DIA) proteomics data, as well as a real-world case study using the UPenn cohort.

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

The `DDA/` directory contains R scripts for reproducing all results from the manuscript's DDA-based proteomics analysis.

### 1. Running Benchmark Methods

- **Script:** `FragPipe.run.r`, `FragPipe.MSstats.r`, and `FragPipe.DEqMS.r` , and `FragPipe.DEP.r`
- **Input Data:** Download the required FragPipe-based protein expression datasets from [Zenodo (DOI:10.5281/zenodo.10953347)](https://zenodo.org/records/10953347).
- **Output Folder:** Create a folder named `FragPipe_results/` to store the method outputs.

> 📌 The same instructions apply to `Maxquant.run.r`, using a `Maxquant_results/` output directory.

### 2. Evaluation of Results

- **Script:** `evaluate.r`
- **Function:** Calculates evaluation metrics across all methods.
- **Input:** Method outputs stored in `FragPipe_results/` or `Maxquant_results/`.
- **Output Folder:** Metrics are saved in `metrics_FragPipe/` (or `metrics_Maxquant/`, which should be created manually).

---

## DIA Analysis

The `DIA/` directory includes scripts to reproduce the DIA-based analyses from Case Studies 1 through 4.

### 1. Running Benchmark Methods

- **Scripts:** 
  - `spectronaut.run.r`, `spectronaut.MSstats.r`, `spectronaut.DEqMS.r`, and `spectronaut.DEP.r`
  - `DIANN.run.r`, `DIANN.MSstats.r`, `DIANN.DEqMS.r`, and `DIANN.DEP.r`
- **Input Data:** Same datasets as above, from [Zenodo (DOI:10.5281/zenodo.10953347)](https://zenodo.org/records/10953347).
- **Output Folders:** 
  - For Spectronaut: `spectronaut_results/`
  - For DIA-NN: `DIANN_results/`

> 📌 Ensure these folders are created before running the scripts.

### 2. Evaluation of Results

- **Script:** `evaluate.r`
- **Function:** Computes performance metrics for each method and saves results.
- **Output Folder:** 
  - For Spectronaut: `metrics_spectronaut/`
  - For DIA-NN: `metrics_DIANN/`

Repeat the same procedure for all four DIA case studies.

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

- Please ensure that all required folders (`*_results/`, `metrics_*`) are created prior to running the scripts.
- All R scripts assume that the required dependencies and packages are installed.

---

## Citation

If you use this repository, please cite the LimROTS manuscript accordingly.

