#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Proteomics Analysis Pipeline
# This script runs the complete analysis pipeline for DIA-NN or Spectronaut data
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
check_and_create_dirs <- function() {
    message("\nChecking and creating required directories:")
    message(paste(rep("-", 50), collapse=""))
    
    # Input directories
    if (!dir.exists(CURRENT_DIRS$input)) {
        stop(sprintf("Input directory not found: %s\nPlease create this directory and add your input files.", CURRENT_DIRS$input))
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
    
    # Verify output directories are writable
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
    
    # Check for required scripts
    required_scripts <- if(ANALYSIS_MODE == "DIANN") {
        c("DIANN.DEP.r", "DIANN.DEqMS.r", "DIANN.MSstats.r")
    } else {
        c("Spectronaut.DEP.r", "Spectronaut.DEqMS.r", "Spectronaut.MSstats.r")
    }
    
    for (script in required_scripts) {
        if (!file.exists(script)) {
            stop(sprintf("Missing required script: %s", script))
        }
        message(sprintf("✓ Found %s", script))
    }
    
    # Check for input data files based on mode
    if (length(list.files(CURRENT_DIRS$input)) == 0) {
        stop(sprintf("No input files found in %s directory", CURRENT_DIRS$input))
    }
    message(sprintf("✓ Input files found in %s", CURRENT_DIRS$input))
    
    message("\nAll required files present!")
}

run_script <- function(script_name, description) {
    message(sprintf("\n%s\n%s\n", description, paste(rep("-", nchar(description)), collapse="")))
    
    if (!file.exists(script_name)) {
        stop(sprintf("Script not found: %s", script_name))
    }
    
    tryCatch({
        source(script_name, echo = TRUE)
        message(sprintf("Completed: %s", description))
    }, error = function(e) {
        message(sprintf("Error in %s: %s", script_name, e$message))
        stop(e)
    })
}

#------------------------------------------------------------------------------
# Package Management
#------------------------------------------------------------------------------
check_required_packages <- function() {
    message("\nChecking for required R packages:")
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
    
    message("\nChecking package installations...")
    missing_packages <- setdiff(required_packages, rownames(installed.packages()))
    
    # Install missing CRAN packages
    if (length(missing_packages) > 0) {
        message("Installing missing packages:", paste(missing_packages, collapse=", "))
        for (pkg in missing_packages) {
            if (!pkg %in% rownames(installed.packages())) {
                # Try BiocManager first
                if (!requireNamespace("BiocManager", quietly = TRUE)) {
                    install.packages("BiocManager", quiet = TRUE)
                }
                tryCatch(
                    BiocManager::install(pkg, update = FALSE, ask = FALSE),
                    error = function(e) {
                        # If BiocManager fails, try CRAN
                        install.packages(pkg, quiet = TRUE)
                    }
                )
            }
        }
    }
    
    message("✓ All required packages are available")
    
    # Load required packages
    message("\nLoading required packages...")
    sapply(required_packages, require, character.only = TRUE, quietly = TRUE)
    message("✓ All packages loaded successfully")
}

#------------------------------------------------------------------------------
# Pipeline Execution Functions
#------------------------------------------------------------------------------
run_analysis <- function(tool_name) {
    script_prefix <- if(ANALYSIS_MODE == "DIANN") "DIANN" else "spectronaut"
    script_name <- sprintf("%s.%s.r", script_prefix, tool_name)
    
    message(sprintf("\nRunning %s analysis...", tool_name))
    tryCatch({
        source(script_name)
        message(sprintf("✓ Completed %s analysis", tool_name))
    }, error = function(e) {
        stop(sprintf("Error in %s analysis: %s", tool_name, e$message))
    })
}

run_complete_pipeline <- function() {
    # Set analysis start time
    start_time <- Sys.time()
    
    message(sprintf("\nStarting %s Analysis Pipeline", ANALYSIS_MODE))
    message(paste(rep("=", 60), collapse=""))
    message(sprintf("Started at: %s", format(start_time, "%Y-%m-%d %H:%M:%S")))
    
    tryCatch({
        # Initialize environment
        check_required_packages()
        check_and_create_dirs()
        check_required_files()
        
        # Display configuration
        message("\nPipeline Configuration:")
        message(sprintf("Analysis Mode: %s", ANALYSIS_MODE))
        message(sprintf("Input Directory: %s", CURRENT_DIRS$input))
        message(sprintf("Results Directory: %s", CURRENT_DIRS$results))
        
        # Execute pipeline steps
        message("\nExecuting pipeline steps:")
        message(paste(rep("-", 60), collapse=""))
        
        # Step 1: Base Statistical Analysis
        message("\n1. Core Statistical Methods")
        if (ANALYSIS_MODE == "DIANN") {
            run_script("DIANN.run.r", "Running DIA-NN analysis (ROTS, LimROTS, t-test, ANOVA, SAM)")
        } else {
            run_script("spectronaut.run.r", "Running Spectronaut analysis (ROTS, LimROTS, t-test, ANOVA, SAM)")
        }
        
        # Step 2: Specialized Tool Analysis
        message("\n2. Specialized Proteomics Tools")
        if (ANALYSIS_MODE == "DIANN") {
            run_script("DIANN.DEP.r", "Running DEP (Differential Enrichment analysis)")
            run_script("DIANN.DEqMS.r", "Running DEqMS (Differential Expression analysis)")
            run_script("DIANN.MSstats.r", "Running MSstats analysis")
        } else {
            run_script("Spectronaut.DEP.r", "Running DEP (Differential Enrichment analysis)")
            run_script("Spectronaut.DEqMS.r", "Running DEqMS (Differential Expression analysis)")
            run_script("Spectronaut.MSstats.r", "Running MSstats analysis")
        }
        
        # Step 3: Evaluation and Visualization
        message("\n3. Results Analysis")
        run_script("evaluate.r", "Running evaluation metrics")
        run_script("plot.r", "Generating visualization plots")
        run_script("sumROC.r", "Generating ROC/PR curves and performance metrics")
        
        # Calculate pipeline duration
        end_time <- Sys.time()
        duration <- difftime(end_time, start_time, units = "mins")
        
        # Print completion summary
        message("\nPipeline Execution Complete")
        message(paste(rep("=", 60), collapse=""))
        message(sprintf("Started:   %s", format(start_time, "%Y-%m-%d %H:%M:%S")))
        message(sprintf("Completed: %s", format(end_time, "%Y-%m-%d %H:%M:%S")))
        message(sprintf("Duration:  %.1f minutes", as.numeric(duration)))
        message(sprintf("Analysis Mode: %s", ANALYSIS_MODE))
        message(sprintf("Output Location: %s", CURRENT_DIRS$results))
        message("\nResults are ready in the following directories:")
        message(sprintf("- Analysis Results: %s", CURRENT_DIRS$results))
        message(sprintf("- Performance Metrics: %s", CURRENT_DIRS$metrics))
        message(sprintf("- Visualization Plots: %s", CURRENT_DIRS$plots))
        message("\nPipeline completed successfully!")
        
    }, error = function(e) {
        message("\nPipeline Error:")
        message(paste(rep("=", 60), collapse=""))
        message(sprintf("Error occurred at: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
        message(sprintf("Error message: %s", e$message))
        message("\nPlease check the error message above and verify:")
        message("1. All required input files are present")
        message("2. Directory permissions are correct")
        message("3. Required R packages are installed")
        stop(e)
    })
}

# Run the complete pipeline
#------------------------------------------------------------------------------
cat("\n")
cat(paste(rep("=", 80), collapse=""))
cat("\n")
cat("LimROTS: Proteomics Analysis Pipeline")
cat("\n")
cat(paste(rep("=", 80), collapse=""))
cat(sprintf("\nStarting %s analysis pipeline...\n", ANALYSIS_MODE))

# Execute the pipeline
run_complete_pipeline()
