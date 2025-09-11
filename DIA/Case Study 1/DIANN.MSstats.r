#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# DIANN.MSstats.r - Differential Expression Analysis using MSstats package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on DIA-NN proteomics data
# using the MSstats package. It processes multiple experiments and performs
# pairwise comparisons between conditions using MSstats-specific methodology.
#
# Required Input Data Structure:
# - DIA-NN output directory containing:
#   * {experiment}/DIANN/report.tsv (DIA-NN output)
#   * {experiment}_DIANN_design.tsv (experimental design)
#   * {experiment}_DIANN_design_msstats.tsv (MSstats annotation)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions should be single characters (e.g., 'A', 'B', 'C')
#
# Output:
# - RDS files containing MSstats analysis results for each comparison
# - Raw MSstats-formatted data saved separately
#------------------------------------------------------------------------------

# Load required libraries
suppressPackageStartupMessages({
    library(stringr)
    library(SummarizedExperiment)
    library(LimROTS)
    library(ROTS)
    library(BiocParallel)
    library(rrcovNA)
    library(imputeLCMD)
    library(limma)
    library(samr)
    library(MSstats)
    library(DEqMS)
    library(matrixStats)
    library(data.table)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
DIANN_DIR <- "../DIANN"  # DIA-NN directory
OUTPUT_DIR <- "DIANN_results/"  # Directory for results
RAW_MSSTATS_DIR <- "rawMSstats/"  # Directory for raw MSstats data

# Create output directories if they don't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)
if (!dir.exists(RAW_MSSTATS_DIR)) dir.create(RAW_MSSTATS_DIR, recursive = TRUE)

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# List all experiments, excluding specific folders
all.exp <- list.files(DIANN_DIR, pattern = "_DIA_DIANN_dlfq.tsv")
all.exp <- str_split_fixed(all.exp, "_DIA_DIANN_", 2)[,1]
all.exp <- unique(all.exp)
#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
    evidence_file <- file.path(DIANN_DIR, paste0(exp_name  , "_DIA"), "DIANN/report.tsv")
    design_file <- file.path(DIANN_DIR, paste0(exp_name, "_DIA_DIANN_design.tsv"))
    msstats_file <- file.path(DIANN_DIR, paste0(exp_name, "_DIA_DIANN_design_msstats.tsv"))
    
    if (!file.exists(evidence_file)) stop("Evidence file not found: ", evidence_file)
    if (!file.exists(design_file)) stop("Design file not found: ", design_file)
    if (!file.exists(msstats_file)) stop("MSstats design file not found: ", msstats_file)
    
    return(list(
        evidence_file = evidence_file,
        design_file = design_file,
        msstats_file = msstats_file
    ))
}


#------------------------------------------------------------------------------
# Main Processing
#------------------------------------------------------------------------------
# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        evidence <- fread(files$evidence_file, sep = "\t", nThread = 15)
        design <- read.csv(files$design_file, sep = "\t")
        annot <- read.csv(files$msstats_file, sep = "\t")
        
        message("  Converting data to MSstats format...")
        # Format data for MSstats
        quant <- DIANNtoMSstatsFormat(
            input = evidence, annot = annot,
            removeFewMeasurements = FALSE, nThread = 15,
            use_log_file = FALSE)
        
        # Save raw MSstats formatted data
        raw_file <- file.path(RAW_MSSTATS_DIR, 
                             paste0(i.exp, "_raw_msstats.rds"))
        saveRDS(quant, raw_file)
        message("  Raw MSstats data saved to: ", raw_file)
        
        # Process data with MSstats
        message("  Running MSstats analysis...")
        processed_data <- dataProcess(quant,
                                    normalization = FALSE, numberOfCores = 15,
                                    use_log_file = FALSE)
    
        # Generate pairwise contrasts
        message("  Generating contrasts...")
        Contrasts <- combn(design$condition, 2, simplify = FALSE)
        Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
        # Remove self-comparisons
        Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
            length(unique(strsplit(x, "")[[1]])) == 1)]
        
        # Process each contrast
        for(i.cnt in Contrasts) {
            tryCatch({
                # Extract conditions
                cnt.cond <- strsplit(i.cnt, "")[[1]]
                message("    Processing contrast: ", cnt.cond[1], " vs ", cnt.cond[2])
                i.Contrasts1 <- str_split(i.cnt, "")[[1]][1]
                i.Contrasts2 <- str_split(i.cnt, "")[[1]][2]
                message("    Processing contrast: ", cnt.cond[1], " vs ", cnt.cond[2])
                
                comparison_r <- levels(processed_data$ProteinLevelData$GROUP)
                comparison_r <- ifelse(comparison_r == i.Contrasts1, -1,
                                       ifelse(comparison_r == i.Contrasts2, 1, 0))
                
                comparison <- matrix(comparison_r, nrow = 1)
                colnames(comparison) <- levels(processed_data$ProteinLevelData$GROUP)
                row.names(comparison) <- paste0(i.Contrasts1, "-", i.Contrasts2)
                
                # Run statistical test
                test.results <- groupComparison(contrast.matrix = comparison,
                                              data = processed_data,
                                              save_fitted_models = FALSE,
                                              numberOfCores = 15,
                                              use_log_file = FALSE)
                
                # Save results
                results_file <- file.path(OUTPUT_DIR,
                                        paste0(i.exp, "_", cnt.cond[1], "vs", 
                                               cnt.cond[2], "_MSstats.rds"))
                saveRDS(test.results$ComparisonResult, results_file)
                message("    Results saved to: ", results_file)
                
                # Clean up memory
                rm(test.results)
                gc()
            }, error = function(e) {
                message("    Error in contrast analysis: ", e$message)
            })
        }
        
        # Clean up memory
        rm(processed_data)
        gc()
        
    }, error = function(e) {
        message("Error in experiment ", i.exp, ": ", e$message)
        stop("Critical error in MSstats analysis - stopping execution")
    })
}

message("\nAnalysis complete!")
      