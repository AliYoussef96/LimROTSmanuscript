
# Multi-method Proteomics Analysis (LimROTS, ROTS, Limma, SAM, t-test, ANOVA, MSstats, DEqMS, DEP)

This folder contains three R script for running multiple statistical analyses on MaxQuant proteomics output. Supported methods include:

`In FragPipe.run.r`

- LimROTS
- ROTS
- Limma
- SAM
- t-test
- ANOVA

`In FragPipe.MSstats.r`

- MSstats

`In FragPipe.DEqMs.r`
  
- DEqMS

`In FragPipe.DEP.r`
  
- DEP

## 📦 Requirements

Install R and the following packages:

```r
install.packages(c("stringr", "imputeLCMD"))
BiocManager::install(c("SummarizedExperiment", "ROTS", "BiocParallel", "LimROTS", "MSstats" , "DEqMS", "DEP" , "limma"))
install.packages("samr", repos = "http://cran.us.r-project.org")
```

To reproduce the results, download the protein expression datasets for the FragPipe software from [https://zenodo.org/records/10953347](https://zenodo.org/records/10953347). You must also create a `FragPipe_results` folder to save the results.

The same procedure applies to `Maxquant.run.r`.

The script `evaluate.r` is used to calculate the evaluation metrics used in the study. It reads the results for each method from the `FragPipe_results` (or `Maxquant_results`) folder, calculates the evaluation metrics, and saves them in a `metrics_FragPipe` folder (which should be created by the user).
