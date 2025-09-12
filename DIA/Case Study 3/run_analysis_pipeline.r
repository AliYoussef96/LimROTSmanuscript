#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Proteomics Analysis Pipeline - Case Study 3
# This script runs the complete analysis pipeline for comparing HE-overlapped 
# fragments (HEof) data from DIA-NN or Spectronaut
#------------------------------------------------------------------------------

# Configuration
#------------------------------------------------------------------------------
ANALYSIS_MODE <- "DIANN"  # Change to "Spectronaut" for Spectronaut analysis

# Define all directory paths
DIRS_CONFIG <- list(
    # Input directories
    input = list(
        diann = "../DIANN",
        spectronaut = "../Spectronaut"
    ),
    
    # Output directories
    output = list(
        diann = list(
            results = "DIANN_results",
            metrics = "metrics_DIANN",
            plots = "Plot_files/DIANN"
        ),
        spectronaut = list(
            results = "Spectronaut_results",
            metrics = "metrics_Spectronaut",
            plots = "Plot_files/Spectronaut"
        )
    )
)

# Set current mode directories
CURRENT_DIRS <- if(ANALYSIS_MODE == "DIANN") {
    list(
        input = DIRS_CONFIG$input$diann,
        results = DIRS_CONFIG$output$diann$results,
        metrics = DIRS_CONFIG$output$diann$metrics,
        plots = DIRS_CONFIG$output$diann$plots
    )
} else {
    list(
        input = DIRS_CONFIG$input$spectronaut,
        results = DIRS_CONFIG$output$spectronaut$results,
        metrics = DIRS_CONFIG$output$spectronaut$metrics,
        plots = DIRS_CONFIG$output$spectronaut$plots
    )
}

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
    
    # Input directory - should point to main DIANN/Spectronaut directories
    if (!dir.exists(CURRENT_DIRS$input)) {
        stop(sprintf("Input directory not found: %s\nPlease ensure the main DIA data directory exists.", CURRENT_DIRS$input))
    }
    message(sprintf("✓ Input directory exists: %s", CURRENT_DIRS$input))
    
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

check_required_files <- function() {
    message("\nChecking for required files:")
    message(paste(rep("-", 50), collapse=""))
    
    # Define required files based on analysis mode
    required_files <- if(ANALYSIS_MODE == "DIANN") {
        c(
            "HEof_n600_DIA_DIANN_dlfq.tsv",
            "HEof_n600_DIA_DIANN_design.tsv",
            "HEof_n600_DIA_DIANN_design_msstats.tsv",
            "HEof_w600_DIA_DIANN_dlfq.tsv", 
            "HEof_w600_DIA_DIANN_design.tsv",
            "HEof_w600_DIA_DIANN_design_msstats.tsv"
        )
    } else {
        c(
            "HEof_n600_DIA_spt_dlfq.tsv",
            "HEof_n600_DIA_spt_design.tsv", 
            "HEof_n600_DIA_spt_design_msstats.tsv",
            "HEof_w600_DIA_spt_dlfq.tsv",
            "HEof_w600_DIA_spt_design.tsv",
            "HEof_w600_DIA_spt_design_msstats.tsv"
        )
    }
    
    # Check required files
    for (file in required_files) {
        file_path <- file.path(CURRENT_DIRS$input, file)
        if (!file.exists(file_path)) {
            stop(sprintf("Missing required file: %s", file_path))
        }
        message(sprintf("✓ Found file: %s", file))
    }
    
    message("\nAll required files found!")
}

run_analysis <- function() {
    message("\nRunning analysis pipeline:")
    message(paste(rep("-", 50), collapse=""))
    
    # Run DEqMS analysis
    message("\n1. Running DEqMS analysis...")
    if(ANALYSIS_MODE == "DIANN") {
        source("HEof_DIANN.r")
    } else {
        source("HEof_Spectronaut.r")
    }
    message("✓ DEqMS analysis completed")
    
    # Run evaluation
    message("\n2. Running evaluation metrics...")
    source("evaluate.r")
    message("✓ Evaluation completed")
    
    # Generate plots
    message("\n3. Generating plots...")
    source("plot.r")
    message("✓ Plots generated")
    
    # Generate ROC analysis
    message("\n4. Generating ROC/PR curves...")
    source("sumROC.r")
    message("✓ ROC analysis completed")
    
    message("\nAnalysis pipeline completed successfully!")
}

#------------------------------------------------------------------------------
# Main Execution
#------------------------------------------------------------------------------
main <- function() {
    # Print header
    message("\nProteomics Analysis Pipeline - Case Study 3")
    message(paste(rep("=", 50), collapse=""))
    message(sprintf("Analysis Mode: %s", ANALYSIS_MODE))
    message("Focus: HE-overlapped fragments (HEof) Analysis")
    
    # Run pipeline steps
    check_and_install_packages()
    check_and_create_dirs()
    check_required_files()
    run_analysis()
    
    # Final message
    message("\nPipeline execution completed!")
    message(sprintf("Results are available in: %s", CURRENT_DIRS$results))
    message(sprintf("Plots are available in: %s", CURRENT_DIRS$plots))
}

# Execute main function
main()
