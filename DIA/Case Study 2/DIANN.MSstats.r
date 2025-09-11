#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# DIANN.MSstats_new.r - Differential Expression Analysis using MSstats package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on DIA-NN proteomics data
# comparing narrow and wide window data using the MSstats package.
# 
#
# Required Input Data Structure:
# - DIANN/{experiment}_n600_DIA/DIANN/report.tsv (Narrow window data)
# - DIANN/{experiment}_w600_DIA/DIANN/report.tsv (Wide window data)
# - DIANN/{experiment}_n600_DIA_DIANN_design_msstats.tsv (Narrow annotations)
# - DIANN/{experiment}_w600_DIA_DIANN_design_msstats.tsv (Wide annotations)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions will be appended with "_Narrow" or "_Wide"
#
# Output:
# - RDS files containing MSstats analysis results comparing window settings
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
# Set paths and experiment configurations
DIANN_DIR <- "../DIANN"
OUTPUT_DIR <- "DIANN_results/"
RAW_MSSTATS_DIR <- "rawMSstats/"

# Create output directories if they don't exist
for (dir in c(OUTPUT_DIR, RAW_MSSTATS_DIR)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# Define experiment names to process
all.exp <- c("HEof")  # Add more experiments as needed

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files for narrow/wide window comparison
validate_input_files <- function(exp_name) {
    # Define file paths for narrow and wide window data
    narrow_evidence <- file.path(DIANN_DIR, paste0(exp_name, "_n600_DIA"), "DIANN/report.tsv")
    wide_evidence <- file.path(DIANN_DIR, paste0(exp_name, "_w600_DIA"), "DIANN/report.tsv")
    narrow_annot <- file.path(DIANN_DIR, paste0(exp_name, "_n600_DIA_DIANN_design_msstats.tsv"))
    wide_annot <- file.path(DIANN_DIR, paste0(exp_name, "_w600_DIA_DIANN_design_msstats.tsv"))
    
    # Check if all required files exist
    if (!file.exists(narrow_evidence)) stop("Narrow evidence file not found: ", narrow_evidence)
    if (!file.exists(wide_evidence)) stop("Wide evidence file not found: ", wide_evidence)
    if (!file.exists(narrow_annot)) stop("Narrow annotation file not found: ", narrow_annot)
    if (!file.exists(wide_annot)) stop("Wide annotation file not found: ", wide_annot)
    
    return(list(
        narrow_evidence = narrow_evidence,
        wide_evidence = wide_evidence,
        narrow_annot = narrow_annot,
        wide_annot = wide_annot
    ))
}

# Function to load and prepare narrow/wide window data
load_window_data <- function(files) {
    message("  Loading narrow window data...")
    narrow_exp <- tryCatch({
        data <- fread(files$narrow_evidence, sep = "\t", nThread = 15)
        # Add window type identifier to condition
        data$R.Condition <- paste0(data$R.Condition, "_Narrow")
        data
    }, error = function(e) {
        stop("Error loading narrow window data: ", e$message)
    })
    
    narrow_annot <- tryCatch({
        data <- read.csv(files$narrow_annot, sep = "\t")
        # Add window type identifier to condition
        data$Condition <- paste0(data$Condition, "_Narrow")
        data
    }, error = function(e) {
        stop("Error loading narrow window annotations: ", e$message)
    })
    
    message("  Loading wide window data...")
    wide_exp <- tryCatch({
        data <- fread(files$wide_evidence, sep = "\t", nThread = 15)
        # Add window type identifier to condition
        data$R.Condition <- paste0(data$R.Condition, "_Wide")
        data
    }, error = function(e) {
        stop("Error loading wide window data: ", e$message)
    })
    
    wide_annot <- tryCatch({
        data <- read.csv(files$wide_annot, sep = "\t")
        # Add window type identifier to condition
        data$Condition <- paste0(data$Condition, "_Wide")
        data
    }, error = function(e) {
        stop("Error loading wide window annotations: ", e$message)
    })
    
    return(list(
        narrow_exp = narrow_exp,
        wide_exp = wide_exp,
        narrow_annot = narrow_annot,
        wide_annot = wide_annot
    ))
}

# Function to generate contrasts for narrow vs wide comparison
generate_contrasts <- function(data) {
    # Extract base conditions (without _Narrow/_Wide suffix)
    base_conditions <- unique(gsub("_(Narrow|Wide)$", "", data$narrow_annot$Condition))
    
    # Generate pairwise contrasts between base conditions
    contrasts <- combn(base_conditions, 2, simplify = FALSE)
    contrasts <- unique(sapply(contrasts, function(x) paste(sort(x), collapse = "")))
    
    # Remove self-comparisons (shouldn't happen but just in case)
    contrasts <- contrasts[!sapply(contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    
    return(contrasts)
}

#------------------------------------------------------------------------------
# Main Processing Loop
#------------------------------------------------------------------------------
message("Starting DIA-NN MSstats analysis for narrow vs wide window comparison...")

for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)

    tryCatch({
        # Validate and load input files
        files <- validate_input_files(i.exp)
        data <- load_window_data(files)
        
        # Combine narrow and wide datasets
        message("  Integrating narrow and wide window datasets...")
        combined_exp <- tryCatch({
            combined <- rbind(data$narrow_exp, data$wide_exp)
            # Use Protein.Ids column as Protein.Names for MSstats compatibility
            combined$Protein.Names <- combined$Protein.Ids
            combined
        }, error = function(e) {
            stop("Error combining expression datasets: ", e$message)
        })
        
        combined_annot <- tryCatch({
            rbind(data$narrow_annot, data$wide_annot)
        }, error = function(e) {
            stop("Error combining annotation datasets: ", e$message)
        })
        
        # Generate contrasts
        message("  Generating pairwise contrasts...")
        contrasts <- generate_contrasts(data)
        message(sprintf("  Found %d contrasts to analyze: %s", 
                       length(contrasts), paste(contrasts, collapse = ", ")))
        
        # Convert DIA-NN format to MSstats format
        message("  Converting data to MSstats format...")
        raw <- tryCatch({
            DIANNtoMSstatsFormat(
                annotation = combined_annot,
                input = combined_exp,
                removeFewMeasurements = FALSE, 
                nThread = 15,
                use_log_file = FALSE
            )
        }, error = function(e) {
            stop("Error converting to MSstats format: ", e$message)
        })
        
        # Save raw MSstats formatted data
        raw_file <- file.path(RAW_MSSTATS_DIR, 
                             paste0("RawMSstats_", i.exp, ".rds"))
        saveRDS(raw, raw_file)
        message("  Raw MSstats data saved to: ", raw_file)
        
        # Process data with MSstats
        message("  Running MSstats data processing...")
        processed_data <- tryCatch({
            dataProcess(raw,
                       normalization = FALSE, 
                       numberOfCores = 20,
                       use_log_file = FALSE)
        }, error = function(e) {
            stop("Error in MSstats data processing: ", e$message)
        })
        
        # Save processed data
        processed_file <- file.path(RAW_MSSTATS_DIR, 
                                   paste0("QuantDataMSstats_", i.exp, ".rds"))
        saveRDS(processed_data, processed_file)
        message("  Processed MSstats data saved to: ", processed_file)
        
        # Clean up memory before contrast analysis
        rm(raw, combined_exp, data)
        gc()
        
        #------------------------------------------------------------------------------
        # Contrast Analysis Loop
        #------------------------------------------------------------------------------
        message("\n  Beginning differential expression analysis...")
        
        for(contrast in contrasts) {
            message(sprintf("\n    Processing contrast: %s", contrast))
            
            tryCatch({
                # Extract individual conditions
                cond1 <- str_split(contrast, "")[[1]][1]
                cond2 <- str_split(contrast, "")[[1]][2]
                
                message(sprintf("    Analyzing: %s vs %s", cond1, cond2))
                
                # Create window-specific condition labels
                cond1_narrow <- paste0(cond1, "_Narrow")
                cond1_wide <- paste0(cond1, "_Wide")
                cond2_narrow <- paste0(cond2, "_Narrow")
                cond2_wide <- paste0(cond2, "_Wide")
                
                # Check if all conditions exist in the data
                available_groups <- levels(processed_data$ProteinLevelData$GROUP)
                required_groups <- c(cond1_narrow, cond1_wide, cond2_narrow, cond2_wide)
                missing_groups <- setdiff(required_groups, available_groups)
                
                if(length(missing_groups) > 0) {
                    message(sprintf("    Skipping contrast: missing groups %s", 
                                   paste(missing_groups, collapse = ", ")))
                    next
                }
                
                # Generate contrast matrix for narrow/wide window comparison
                message("    Setting up comparison matrix...")
                comparison_r <- rep(0, length(available_groups))
                names(comparison_r) <- available_groups
                
                # Set up contrast: (cond1_narrow + cond1_wide) vs (cond2_narrow + cond2_wide)
                comparison_r[cond1_narrow] <- 0.5
                comparison_r[cond1_wide] <- 0.5
                comparison_r[cond2_narrow] <- -0.5
                comparison_r[cond2_wide] <- -0.5
                
                comparison <- matrix(comparison_r, nrow = 1)
                colnames(comparison) <- available_groups
                row.names(comparison) <- paste0(cond1, "-", cond2)
                
                # Perform group comparison
                message("    Running MSstats group comparison...")
                test_results <- tryCatch({
                    groupComparison(contrast.matrix = comparison,
                                   data = processed_data,
                                   numberOfCores = 20,
                                   save_fitted_models = FALSE,
                                   use_log_file = FALSE)
                }, error = function(e) {
                    stop("Error in group comparison: ", e$message)
                })
                
                # Save results
                results_file <- file.path(OUTPUT_DIR,
                                         sprintf("MSstats_%s_%s.rds", i.exp, contrast))
                saveRDS(test_results$ComparisonResult, results_file)
                message(sprintf("    Results saved to: %s", results_file))
                
                # Report summary statistics
                if(!is.null(test_results$ComparisonResult)) {
                    sig_proteins <- sum(test_results$ComparisonResult$adj.pvalue < 0.05, na.rm = TRUE)
                    total_proteins <- nrow(test_results$ComparisonResult)
                    message(sprintf("    Found %d/%d significantly different proteins (FDR < 0.05)", 
                                   sig_proteins, total_proteins))
                }
                
                # Clean up memory
                rm(test_results)
                gc()
                
            }, error = function(e) {
                message(sprintf("    Error processing contrast %s: %s", contrast, e$message))
            })
        }
        
        # Clean up memory before next experiment
        rm(processed_data, combined_annot)
        gc()
        
        message(sprintf("\nCompleted processing experiment: %s", i.exp))
        
    }, error = function(e) {
        message(sprintf("Error in experiment %s: %s", i.exp, e$message))
    })
}

#------------------------------------------------------------------------------
# Final Summary
#------------------------------------------------------------------------------
message("DIA-NN MSstats analysis completed!")
message("Results stored in: ", normalizePath(OUTPUT_DIR))
message("Raw data stored in: ", normalizePath(RAW_MSSTATS_DIR))
message("Analysis focused on narrow vs wide window comparison using MSstats methodology")
