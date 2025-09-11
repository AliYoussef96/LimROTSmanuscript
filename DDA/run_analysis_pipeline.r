#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Proteomics Analysis Pipeline
# This script runs the complete analysis pipeline for MaxQuant or FragPipe data
#------------------------------------------------------------------------------

# Configuration
#------------------------------------------------------------------------------
ANALYSIS_MODE <- "FragPipe"  # Change to "MaxQuant" for MaxQuant analysis

# Define all directory paths
DIRS_CONFIG <- list(
    # Input directories
    input = list(
        maxquant = "../Maxquant",
        fragpipe = "../FragPipe"
    ),
    
    # Output directories
    output = list(
        maxquant = list(
            results = "Maxquant_results",
            metrics = "metrics_Maxquant",
            plots = "Plot_files/MaxQuant"
        ),
        fragpipe = list(
            results = "FragPipe_results",
            metrics = "metrics_FragPipe",
            plots = "Plot_files/FragPipe"
        )
    )
)

# Set current mode directories
CURRENT_DIRS <- if(ANALYSIS_MODE == "MaxQuant") {
    list(
        input = DIRS_CONFIG$input$maxquant,
        results = DIRS_CONFIG$output$maxquant$results,
        metrics = DIRS_CONFIG$output$maxquant$metrics,
        plots = DIRS_CONFIG$output$maxquant$plots
    )
} else {
    list(
        input = DIRS_CONFIG$input$fragpipe,
        results = DIRS_CONFIG$output$fragpipe$results,
        metrics = DIRS_CONFIG$output$fragpipe$metrics,
        plots = DIRS_CONFIG$output$fragpipe$plots
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

check_required_packages <- function() {
    required_packages <- c(
        "stringr", "SummarizedExperiment", "LimROTS", "ROTS",
        "BiocParallel", "rrcovNA", "imputeLCMD", "limma", "samr",
        "DEP", "DEqMS", "MSstats", "dplyr", "ggplot2", 
        "ggstatsplot", "ggsignif", "pROC", "PRROC", "gridExtra", "ggsci"
    )
    
    missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]
    
    if (length(missing_packages) > 0) {
        stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
    }
}

#------------------------------------------------------------------------------
# Main Pipeline
#------------------------------------------------------------------------------
main <- function() {
    # Print header with timestamp
    timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
    message(sprintf("\nProteomics Analysis Pipeline (%s mode)", ANALYSIS_MODE))
    message(sprintf("Started: %s", timestamp))
    message(paste(rep("=", 60), collapse=""))
    
    # Initialize pipeline
    message("\nInitializing pipeline...")
    tryCatch({
        # Check package requirements
        message("\nChecking R package requirements...")
        check_required_packages()
        
        # Setup directories
        check_and_create_dirs()
        
        # Verify input files
        message("\nVerifying input files...")
        check_required_files()
        
        # Display pipeline configuration
        message("\nPipeline Configuration:")
        message(sprintf("- Analysis Mode: %s", ANALYSIS_MODE))
        message(sprintf("- Input Directory: %s", CURRENT_DIRS$input))
        message(sprintf("- Results Directory: %s", CURRENT_DIRS$results))
        message(sprintf("- Metrics Directory: %s", CURRENT_DIRS$metrics))
        message(sprintf("- Plots Directory: %s", CURRENT_DIRS$plots))
        
        # Execute pipeline steps
        message("\nExecuting pipeline steps:")
        message(paste(rep("-", 60), collapse=""))
        
        # Step 1: Base Statistical Analysis
        message("\n1. Core Statistical Methods")
        if (ANALYSIS_MODE == "MaxQuant") {
            run_script("Maxquant.run.r", "Running MaxQuant analysis (ROTS, LimROTS, t-test, ANOVA, SAM)")
        } else {
            run_script("FragPipe.run.r", "Running FragPipe analysis (ROTS, LimROTS, t-test, ANOVA, SAM)")
        }
        
        # Step 2: Specialized Tool Analysis
        message("\n2. Specialized Proteomics Tools")
        if (ANALYSIS_MODE == "MaxQuant") {
            run_script("Maxquant.DEP.r", "Running DEP (Differential Enrichment analysis)")
            run_script("Maxquant.DEqMS.r", "Running DEqMS (Differential Expression analysis)")
            run_script("Maxquant.MSstats.r", "Running MSstats analysis")
        } else {
            run_script("FragPipe.DEP.r", "Running DEP (Differential Enrichment analysis)")
            run_script("FragPipe.DEqMS.r", "Running DEqMS (Differential Expression analysis)")
            run_script("FragPipe.MSstats.r", "Running MSstats analysis")
        }
        
        # Step 3: Evaluation and Visualization
        message("\n3. Results Analysis")
        run_script("evaluate.r", "Running evaluation metrics")
        run_script("sumROC.r", "Generating ROC/PR curve analysis") 
        run_script("plot.r", "Generating visualization plots")
    
        # Print completion summary
        end_timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
        message("\nPipeline Execution Summary")
        message(paste(rep("=", 60), collapse=""))
        message(sprintf("Started: %s", timestamp))
        message(sprintf("Completed: %s", end_timestamp))
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

# Run the pipeline
tryCatch({
    main()
}, error = function(e) {
    message("\nFatal error occurred. Pipeline execution stopped.")
    quit(status = 1)
})
