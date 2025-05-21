
# Multi-method Proteomics Analysis (LimROTS, ROTS, Limma, SAM, t-test, ANOVA, MSstats, DEqMS, DEP)

This folder contains three R script for running multiple statistical analyses for each Spectronaut and DIANN proteomics output (case studies 1, 2, 3, and 4). Supported methods include:

`In DIANN.run.r` and `In Spectronaut.run.r`

- LimROTS
- ROTS
- Limma
- SAM
- t-test
- ANOVA

`In DIANN.MSstats.r` and `In Spectronaut.MSstats.r`

- MSstats

`In DIANN.DEqMs.r` and `In Spectronaut.DEqMs.r`
  
- DEqMS

`In DIANN.DEP.r` and `In Spectronaut.DEP.r`
  
- DEP

## Requirements

Install R and the following packages:

```r
install.packages(c("stringr", "imputeLCMD"))
BiocManager::install(c("SummarizedExperiment", "ROTS", "BiocParallel", "LimROTS", "MSstats" , "DEqMS", "DEP" , "limma"))
install.packages("samr", repos = "http://cran.us.r-project.org")
```

## To reproduce the results

Download the protein expression datasets for the DIANN software from [https://zenodo.org/records/10953347](https://zenodo.org/records/10953347). You must also create a `DIANN_results` folder to save the results.

For `MSstats` and `DEqMS` you should also download the raw files from: [https://zenodo.org/records/10482353](https://zenodo.org/records/10482353).

The same procedure applies to `Spectronaut.run.r`.

The script `evaluate.r` is used to calculate the evaluation metrics used in the study. It reads the results for each method from the `DIANN_results` (or `Spectronaut_results`) folder, calculates the evaluation metrics, and saves them in a `metrics_DIANN` folder (which should be created by the user).
