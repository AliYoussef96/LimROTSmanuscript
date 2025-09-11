#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Proteomics Analysis Pipeline - Case Study 4
# This script runs the complete analysis pipeline for differential expression
# analysis using HEof approach
#------------------------------------------------------------------------------

# Configuration
#------------------------------------------------------------------------------
ANALYSIS_MODE <- "HEof"  # Analysis mode for Case Study 4

# Define output directory paths
DIRS_CONFIG <- list(
    output = list(
        results = "results",
        metrics = "metrics",
        plots = "Plot_files"
    )
)

# Set current directories
CURRENT_DIRS <- list(
    results = DIRS_CONFIG$output$results,
    metrics = DIRS_CONFIG$output$metrics,
    plots = DIRS_CONFIG$output$plots
)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
check_and_install_packages <- function() {
    message("\nChecking and installing required packages:")
    message(paste(rep("-", 50), collapse=""))
    
    required_packages <- c(
        # Core analysis packages
        "stringr", "SummarizedExperiment", "LimROTS", "ROTS",
        "BiocParallel", "rrcovNA", "imputeLCMD", "limma", "samr",
        "DEP", "DEqMS", "MSstats",
        
        # Data manipulation and visualization
        "dplyr", "tidyr", "ggplot2", "patchwork",
        "ggrepel", "pheatmap", "RColorBrewer",
        
        # ROC analysis packages
        "pROC", "PRROC", "gridExtra", "ggsci", "caret"
    )
    
    missing_packages <- setdiff(required_packages, rownames(installed.packages()))
    
    if(length(missing_packages) > 0) {
        message("Installing missing packages...")
        
        # Separate CRAN and Bioconductor packages
        bioc_packages <- c("SummarizedExperiment", "LimROTS", "ROTS", "BiocParallel", 
                          "rrcovNA", "imputeLCMD", "limma", "samr", "DEP", "DEqMS", "MSstats")
        cran_packages <- setdiff(missing_packages, bioc_packages)
        missing_bioc <- intersect(missing_packages, bioc_packages)
        
        # Install CRAN packages
        if(length(cran_packages) > 0) {
            install.packages(cran_packages, quiet = TRUE)
        }
        
        # Install Bioconductor packages
        if(length(missing_bioc) > 0) {
            if(!requireNamespace("BiocManager", quietly = TRUE)) {
                install.packages("BiocManager", quiet = TRUE)
            }
            BiocManager::install(missing_bioc, quiet = TRUE)
        }
    }
    
    # Load all packages
    suppressMessages({
        sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
    })
    
    message("✓ All packages loaded successfully!")
}

check_and_create_dirs <- function() {
    message("\nChecking and creating required directories:")
    message(paste(rep("-", 50), collapse=""))
    
    message("✓ No input directories needed (using built-in UPS1.Case4 dataset)")
    
    # Create output directories if they don't exist
    output_dirs <- c(CURRENT_DIRS$results, CURRENT_DIRS$metrics, CURRENT_DIRS$plots)
    for (dir in output_dirs) {
        if (!dir.exists(dir)) {
            dir.create(dir, recursive = TRUE)
            message(sprintf("Created output directory: %s", dir))
        } else {
            message(sprintf("✓ Output directory exists: %s", dir))
        }
    }
    
    # Verify all output directories are writable
    for (dir in output_dirs) {
        if (!file.access(dir, mode=2) == 0) {
            stop(sprintf("Directory not writable: %s\nPlease check permissions.", dir))
        }
    }
    
    message("\nAll directory checks completed successfully!")
}

check_data_availability <- function() {
    message("\nChecking data availability:")
    message(paste(rep("-", 50), collapse=""))
    
    # Try to load UPS1.Case4 dataset
    tryCatch({
        # Load LimROTS package which should contain the dataset
        if(!require("LimROTS", quietly = TRUE)) {
            stop("LimROTS package not available - required for UPS1.Case4 dataset")
        }
        
        data(UPS1.Case4)
        if (!exists("UPS1.Case4")) {
            stop("UPS1.Case4 dataset not found")
        }
        message("✓ Successfully loaded UPS1.Case4 dataset")
        
        # Validate dataset structure
        if (!all(c("tool", "Conc.", "fake.batch") %in% colnames(colData(UPS1.Case4)))) {
            stop("Missing required metadata columns in UPS1.Case4")
        }
        message("✓ Dataset structure validated")
        
    }, error = function(e) {
        stop("Failed to load UPS1.Case4 dataset: ", e$message)
    })
    
    message("\nData check completed successfully!")
}

run_analysis <- function() {
    message("\nRunning analysis pipeline:")
    message(paste(rep("-", 50), collapse=""))
    
    # Run HEof analysis
    message("\n1. Running HEof analysis...")
    source("HEof.r")
    message("✓ HEof analysis completed")
    
    # Run evaluation
    message("\n2. Running evaluation metrics...")
    source("evaluate.r")
    message("✓ Evaluation completed")
    
    # Generate plots
    message("\n3. Generating plots...")
    source("plot.r")
    message("✓ Plots generated")
    
    # Note: ROC analysis may not be applicable for this case study
    # as it focuses on HEof approach rather than method comparison
    # If sumROC.r exists, uncomment the following lines:
    # message("\n4. Generating ROC/PR curves...")
    # source("sumROC.r")
    # message("✓ ROC analysis completed")
    
    message("\nAnalysis pipeline completed successfully!")
}

#------------------------------------------------------------------------------
# Main Execution
#------------------------------------------------------------------------------
main <- function() {
    # Print header
    message("\nProteomics Analysis Pipeline - Case Study 4")
    message(paste(rep("=", 50), collapse=""))
    message(sprintf("Analysis Mode: %s", ANALYSIS_MODE))
    message("Focus: HE-overlapped fragments (HEof) Analysis")
    
    # Run pipeline steps
    check_and_install_packages()
    check_and_create_dirs()
    check_data_availability()
    run_analysis()
    
    # Final message
    message("\nPipeline execution completed!")
    message(sprintf("Results are available in: %s", CURRENT_DIRS$results))
    message(sprintf("Plots are available in: %s", CURRENT_DIRS$plots))
}

# Execute main function
main()
