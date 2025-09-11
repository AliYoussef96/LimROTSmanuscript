#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# DIANN.DEP.r - Differential Expression Analysis using DEP package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on DIA-NN proteomics data
# using the DEP package. It processes multiple experiments and performs pairwise
# comparisons between conditions.
#
# Required Input Data Structure:
# - DIA-NN output directory containing:
#   * {experiment}_DIA_DIANN_dlfq.tsv (intensity data)
#   * {experiment}_LFQ_DIANN_design.tsv (experimental design)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions should be single characters (e.g., 'A', 'B', 'C')
#
# Output:
# - RDS files containing differential expression results for each comparison
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
    library(multiUS)
    library(DEP)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
INPUT_DIR <- "../DIANN/"  # Directory containing DIA-NN output
OUTPUT_DIR <- "DIANN_results/"  # Directory for results
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
  intensity_file <- file.path(INPUT_DIR, paste0(exp_name, "_DIA_DIANN_dlfq.tsv"))
  design_file <- file.path(INPUT_DIR, paste0(exp_name, "_DIA_DIANN_design.tsv"))
  
  if (!file.exists(intensity_file)) stop("Intensity file not found: ", intensity_file)
  if (!file.exists(design_file)) stop("Design file not found: ", design_file)
  
  return(list(intensity_file = intensity_file, design_file = design_file))
}

# List all files in the directory and extract experiment names
all.exp <- list.files(INPUT_DIR, pattern = '_DIA_DIANN_dlfq')
all.exp <- unique(str_split_fixed(all.exp, "_DIA_DIANN_", 2)[,1])

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------

# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        exp.mar <- read.csv(files$intensity_file, sep = "\t")
        design <- read.csv(files$design_file, sep = "\t")
    }, error = function(e) {
        message("Error processing experiment ", i.exp, ": ", e$message)
        stop("Critical error in DEP analysis - stopping execution")
    })
    
    # Generate pairwise contrasts
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    
    # Format intensity data
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    
    # Clean sample names by removing underscores
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Ensure data and design matrix match
    exp.mar <- exp.mar[, colnames(exp.mar) %in% design$sample_name]
    design <- design[design$sample_name %in% colnames(exp.mar), ]
    exp.mar <- exp.mar[, design$sample_name]
    
    # Process each contrast
    for(i.Contrasts in Contrasts) {
        message("  Processing contrast: ", i.Contrasts)
        
        # Extract conditions
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        # Count replicates
        select.col.Contrasts1 <- sum(design$condition == i.Contrasts1)
        select.col.Contrasts2 <- sum(design$condition == i.Contrasts2)
        
        # Only process if both conditions have sufficient replicates
        if(select.col.Contrasts1 > 1 && select.col.Contrasts2 > 1) {
            # Subset data for current contrast
            design.temp <- design[design$condition %in% c(i.Contrasts1, i.Contrasts2), ]
            df.temp <- exp.mar[, colnames(exp.mar) %in% design.temp$sample_name]
            df.temp <- df.temp[, design.temp$sample_name]
            
            # Data preprocessing
            # Convert zeros to NA and log2 transform
            df.temp[df.temp == 0] <- NA
            df.temp <- log2(df.temp + 1)
            # Impute missing values using sequential KNN
            df.temp <- impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = FALSE, check.rows = FALSE)
            
            # Prepare experimental design
            experimental_design <- design.temp[, c(3,4,5)]
            colnames(experimental_design)[1] <- "label"
            df.temp$name <- row.names(df.temp)
            df.temp$ID <- row.names(df.temp)
            experimental_design$replicate <- seq(1, nrow(experimental_design))
            
            # Create SummarizedExperiment object
            data_se <- make_se(df.temp, seq(1, ncol(df.temp)-2), experimental_design)
            
            # Clean condition names
            data_se$label <- make.names(data_se$label)
            data_se$condition <- make.names(data_se$condition)
            
            # Perform differential expression analysis
            data_diff_manual <- test_diff(data_se, 
                                        control = i.Contrasts1,
                                        test = i.Contrasts2,
                                        design_formula = formula("~ 0 + condition"))
            
            # Extract and format results
            data_diff_manual.df <- data.frame(elementMetadata(data_diff_manual))
            data_diff_manual.df <- data_diff_manual.df[, c(5,6,7)]
            colnames(data_diff_manual.df) <- str_remove(
                colnames(data_diff_manual.df),
                paste0(i.Contrasts2, "_vs_", i.Contrasts1)
            )
            colnames(data_diff_manual.df) <- str_remove(
                colnames(data_diff_manual.df),
                fixed("_")
            )
            
            # Save results
            output_file <- file.path(OUTPUT_DIR, 
                                   paste0("DEP_", i.exp, "_", i.Contrasts, ".rds"))
            saveRDS(data_diff_manual.df, output_file)
            message("    Results saved to: ", output_file)
        } else {
            message("    Skipping contrast due to insufficient replicates")
        }
    }
}


message("Analysis complete!")
