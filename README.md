# LimROTSmanuscript

## DDA Analysis

Inside the `DDA` folder, R scripts are available to reproduce all the results presented in the manuscript.

First, `FragPipe.run.r` contains the necessary code to run all the methods benchmarked in the study. To reproduce the results, download the protein expression datasets for the FragPipe software from [https://zenodo.org/records/10953347](https://zenodo.org/records/10953347). You must also create a `FragPipe_results` folder to save the results.

The same procedure applies to `Maxquant.run.r`.

The script `evaluate.r` is used to calculate the evaluation metrics used in the study. It reads the results for each method from the `FragPipe_results` (or `Maxquant_results`) folder, calculates the evaluation metrics, and saves them in a `metrics_FragPipe` folder (which should be created by the user).

## DIA Analysis

Inside the `DIA` folder, R scripts are available to reproduce all the results for Case Studies 1 to 4 presented in the manuscript.

Similar to the DDA analysis, `spectronaut.run.r` (and `DIANN.run.r`) contain the necessary code to run all the methods benchmarked in the study. To reproduce the results, download the protein expression datasets for the FragPipe software from [https://zenodo.org/records/10953347](https://zenodo.org/records/10953347). You must also create a `spectronaut_results` folder to save the results. The same procedure applies to `DIANN.run.r`.

The script `evaluate.r` is used to calculate the evaluation metrics used in the study. It reads the results for each method from the `spectronaut_results` (or `DIANN_results`) folder, calculates the evaluation metrics, and saves them in a `metrics_spectronaut` folder (which should be created by the user).

The same procedure applies to Case Studies 2, 3, and 4.

