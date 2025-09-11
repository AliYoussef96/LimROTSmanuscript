#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# FragPipe.MSstats.r - Differential Expression Analysis using MSstats package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on FragPipe proteomics data
# using the MSstats package. It processes multiple experiments and performs
# pairwise comparisons between conditions using MSstats-specific methodology.
#
# Required Input Data Structure:
# - FragPipe output directory containing:
#   * {experiment}/FragPipe/TOP0/combined_protein.tsv (protein quantification)
#   * {experiment}/FragPipe/TOP0/MSstats.csv (MSstats evidence file)
#   * {experiment}_FragPipe_design.tsv (experimental design)
#   * {experiment}_FragPipe_design_msstats.tsv (MSstats annotation)
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
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
INPUT_DIR <- "FragPipe"  # Base directory for experiments
OUTPUT_DIR <- "FragPipe_results/"  # Directory for results
RAW_MSSTATS_DIR <- "rawMSstats/"  # Directory for raw MSstats data

# Create output directories if they don't exist
for(dir in c(OUTPUT_DIR, RAW_MSSTATS_DIR)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
  protein_file <- file.path(INPUT_DIR,  paste0(exp_name , "_LFQ") , "FragPipe/TOP0/combined_protein.tsv")
  msstats_file <-  file.path(INPUT_DIR,  paste0(exp_name , "_LFQ") , "FragPipe/TOP0/MSstats.csv")
  design_file <- file.path(INPUT_DIR, paste0(exp_name, "_LFQ_FragPipe_design.tsv"))
  msstats_design_file <- file.path(INPUT_DIR, 
                                    paste0(exp_name, "_LFQ_FragPipe_design_msstats.tsv"))
    
    for(file in list(protein_file, msstats_file, design_file, msstats_design_file)) {
        if (!file.exists(file)) stop("Required file not found: ", file)
    }
    
    return(list(
        protein_file = protein_file,
        msstats_file = msstats_file,
        design_file = design_file,
        msstats_design_file = msstats_design_file
    ))
}

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# List experiment folders (excluding specific ones)
all.exp <- list.files(INPUT_DIR,  pattern = '_LFQ_FragPipe_dlfq_pro_intensity')
all.exp <- str_split_fixed(all.exp, "_LFQ_FragPipe_", 2)[,1]
all.exp <- unique(all.exp)


# Set random seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------
# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        
        # Read all required files
        exp.mar <- read.csv(files$protein_file, sep = "\t")
        design <- read.csv(files$design_file, sep = "\t")
        evidence <- read.csv(files$msstats_file)
        annot <- read.csv(files$msstats_design_file, sep = "\t")
    }, error = function(e) {
        message("Error reading files for experiment ", i.exp, ": ", e$message)
        stop("Critical error in file reading - stopping execution")
    })
    
    # Generate all valid pairwise condition contrasts
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    
    message("  Processing ", length(Contrasts), " pairwise comparisons")
    
    # Format protein data
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    
    # Clean sample names
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Process data with MSstats
    message("  Converting FragPipe data to MSstats format")
    raw <- tryCatch({
        FragPipetoMSstatsFormat(
            annotation = annot, 
            input = evidence, nThread = 10, removeFewMeasurements = FALSE,
            use_log_file = FALSE
        )
    }, error = function(e) {
        message("    Error in MSstats format conversion: ", e$message)
        stop("Critical error in MSstats format conversion - stopping execution")
    })
    
    if (!is.null(raw)) {
        # Save raw MSstats-formatted data
        raw_file <- file.path(RAW_MSSTATS_DIR, paste0("RawMSstats_", i.exp, ".rds"))
        saveRDS(raw, raw_file)
        message("    Raw MSstats data saved to: ", raw_file)
    
        # Preprocess data for differential testing
        message("  Preprocessing quantitative data")
        QuantData <- tryCatch({
            dataProcess(raw,  numberOfCores = 15, 
                        normalization = FALSE)
        }, error = function(e) {
            message("    Error in data preprocessing: ", e$message)
            stop("Critical error in data preprocessing - stopping execution")
        })
        
        if (!is.null(QuantData)) {
            # Process each contrast
            for(i.Contrasts in Contrasts) {
                message("    Analyzing contrast: ", i.Contrasts)
                
                # Extract contrast conditions
                i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
                i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
                
                # Build contrast matrix
                tryCatch({
                    comparison_r <- levels(QuantData$ProteinLevelData$GROUP)
                    comparison_r <- ifelse(comparison_r == i.Contrasts1, -1,
                                         ifelse(comparison_r == i.Contrasts2, 1, 0))
                    
                    comparison <- matrix(comparison_r, nrow = 1)
                    colnames(comparison) <- levels(QuantData$ProteinLevelData$GROUP)
                    row.names(comparison) <- paste0(i.Contrasts1, "-", i.Contrasts2)
                    
                    # Run group comparison with MSstats
                    message("      Running MSstats group comparison")
                    testResultOneComparison <- groupComparison(
                        contrast.matrix = comparison, 
                        data = QuantData, numberOfCores = 15, use_log_file = FALSE
                    )
                    
                    # Save results
                    output_file <- file.path(OUTPUT_DIR, 
                                           paste0("MSstats_", i.exp, "_", 
                                                 i.Contrasts, ".rds"))
                    saveRDS(testResultOneComparison$ComparisonResult, output_file)
                    message("      Results saved to: ", output_file)
                    
                    # Clean up memory
                    rm(testResultOneComparison)
                    gc()
                    
                }, error = function(e) {
                    message("      Error in contrast analysis: ", e$message)
                    stop("Critical error in MSstats contrast analysis - stopping execution")
                })
            }
            
            # Clean up memory
            rm(raw, evidence, QuantData)
            gc()
        }
    }
}

message("Analysis complete!")
