#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Proteomics Analysis Pipeline - Case Study 2
# This script runs the complete analysis pipeline for comparing narrow and wide
# window DIA data from DIA-NN or Spectronaut
#------------------------------------------------------------------------------

# Configuration
#------------------------------------------------------------------------------
ANALYSIS_MODE <- "DIANN"  # Change to "Spectronaut" for Spectronaut analysis

# Define all directory paths
DIRS_CONFIG <- list(
    # Input directories
    input = list(
        diann = list(
            narrow = "../DIANN/HEof_n600_DIA/DIANN",
            wide = "../DIANN/HEof_w600_DIA/DIANN",
            config = "../DIANN"
        ),
        spectronaut = list(
            narrow = "../Spectronaut/HEof_n600_DIA/Spectronaut",
            wide = "../Spectronaut/HEof_w600_DIA/Spectronaut",
            config = "../Spectronaut"
        )
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
    
    # Input directories
    for (input_dir in c(CURRENT_DIRS$input$narrow, CURRENT_DIRS$input$wide, CURRENT_DIRS$input$config)) {
        if (!dir.exists(input_dir)) {
            stop(sprintf("Input directory not found: %s\nPlease create this directory and add your input files.", input_dir))
        }
        message(sprintf("✓ Input directory exists: %s", input_dir))
    }
    
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
        list(
            narrow = c("report.tsv"),
            wide = c("report.tsv"),
            config = c(
                "HEof_n600_DIA_DIANN_design.tsv",
                "HEof_n600_DIA_DIANN_design_msstats.tsv",
                "HEof_w600_DIA_DIANN_design.tsv",
                "HEof_w600_DIA_DIANN_design_msstats.tsv"
            )
        )
    } else {
        list(
            narrow = c("Report.tsv"),
            wide = c("Report.tsv"),
            config = c(
                "HEof_n600_DIA_spt_design.tsv",
                "HEof_n600_DIA_spt_design_msstats.tsv",
                "HEof_w600_DIA_spt_design.tsv",
                "HEof_w600_DIA_spt_design_msstats.tsv"
            )
        )
    }
    
    # Check narrow window files
    for (file in required_files$narrow) {
        file_path <- file.path(CURRENT_DIRS$input$narrow, file)
        if (!file.exists(file_path)) {
            stop(sprintf("Missing required narrow window file: %s", file_path))
        }
        message(sprintf("✓ Found narrow window file: %s", file))
    }
    
    # Check wide window files
    for (file in required_files$wide) {
        file_path <- file.path(CURRENT_DIRS$input$wide, file)
        if (!file.exists(file_path)) {
            stop(sprintf("Missing required wide window file: %s", file_path))
        }
        message(sprintf("✓ Found wide window file: %s", file))
    }
    
    # Check config files (design files in main directory)
    for (file in required_files$config) {
        file_path <- file.path(CURRENT_DIRS$input$config, file)
        if (!file.exists(file_path)) {
            stop(sprintf("Missing required config file: %s", file_path))
        }
        message(sprintf("✓ Found config file: %s", file))
    }
    
    message("\nAll required files found!")
}

run_analysis <- function() {
    message("\nRunning analysis pipeline:")
    message(paste(rep("-", 50), collapse=""))
    
    # Run DEqMS analysis
    message("\n1. Running DEqMS analysis...")
    if(ANALYSIS_MODE == "DIANN") {
        source("DIANN.DEqMS.r")
    } else {
        source("Spectronaut.DEqMS.r")
    }
    message("✓ DEqMS analysis completed")
    
    # Run MSstats analysis
    message("\n2. Running MSstats analysis...")
    if(ANALYSIS_MODE == "DIANN") {
        source("DIANN.MSstats.r")
    } else {
        source("Spectronaut.MSstats.r")
    }
    message("✓ MSstats analysis completed")
    
    # Run evaluation
    message("\n3. Running evaluation metrics...")
    source("evaluate.r")
    message("✓ Evaluation completed")
    
    # Generate plots
    message("\n4. Generating plots...")
    source("plot.r")
    message("✓ Plots generated")
    
    # Generate ROC analysis
    message("\n5. Generating ROC/PR curves...")
    source("sumROC.r")
    message("✓ ROC analysis completed")
    
    message("\nAnalysis pipeline completed successfully!")
}

#------------------------------------------------------------------------------
# Main Execution
#------------------------------------------------------------------------------
main <- function() {
    # Print header
    message("\nProteomics Analysis Pipeline - Case Study 2")
    message(paste(rep("=", 50), collapse=""))
    message(sprintf("Analysis Mode: %s", ANALYSIS_MODE))
    
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
