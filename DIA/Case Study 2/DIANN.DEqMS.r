#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# DIANN.DEqMS.r - Differential Expression Analysis using DEqMS package
#------------------------------------------------------------------------------
# This script performs differential expression analysis comparing narrow and wide
# window DIA experiments using the DEqMS package.
#
# Required Input Data Structure:
# - DIA-NN output directories:
#   * HEof_n600_DIA/DIANN/report.tsv (Narrow window data)
#   * HEof_w600_DIA/DIANN/report.tsv (Wide window data)
#   * HEof_n600_DIA_DIANN_dlfq.tsv (Narrow window quantification)
#   * HEof_w600_DIA_DIANN_dlfq.tsv (Wide window quantification)
#   * HEof_n600_DIA_DIANN_design.tsv (Narrow window design)
#   * HEof_w600_DIA_DIANN_design.tsv (Wide window design)
#   * DIANN/*_DIANN_design_msstats.tsv (MSstats annotations)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions will be appended with "_Narrow" or "_Wide"
#
# Output:
# - RDS files containing DEqMS analysis results comparing window settings
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
    library(multiUS)
    library(dplyr)
    library(data.table)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
DIANN_DIR <- "../DIANN"  # Base directory for DIA-NN data
OUTPUT_DIR <- "DIANN_results/"  # Directory for results

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# Define experiment pairs (narrow and wide window experiments)
all.exp.pairs <- list(
    "HEof" = list(narrow = "HEof_n600", wide = "HEof_w600")
)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files for narrow/wide window pair
validate_input_files <- function(exp_narrow, exp_wide) {
    # Narrow window files
    peptide_file_narrow <- file.path(DIANN_DIR, paste0(exp_narrow, "_DIA"), "DIANN/report.tsv")
    design_file_narrow <- file.path(DIANN_DIR, paste0(exp_narrow, "_DIA_DIANN_design.tsv"))
    exp_file_narrow <- file.path(DIANN_DIR, paste0(exp_narrow, "_DIA_DIANN_dlfq.tsv"))
    msstats_design_narrow <- file.path(DIANN_DIR, paste0(exp_narrow, "_DIA_DIANN_design_msstats.tsv"))
    
    # Wide window files
    peptide_file_wide <- file.path(DIANN_DIR, paste0(exp_wide, "_DIA"), "DIANN/report.tsv")
    design_file_wide <- file.path(DIANN_DIR, paste0(exp_wide, "_DIA_DIANN_design.tsv"))
    exp_file_wide <- file.path(DIANN_DIR, paste0(exp_wide, "_DIA_DIANN_dlfq.tsv"))
    msstats_design_wide <- file.path(DIANN_DIR, paste0(exp_wide, "_DIA_DIANN_design_msstats.tsv"))
    
    # Check all files exist
    files_to_check <- list(
        peptide_file_narrow, design_file_narrow, exp_file_narrow, msstats_design_narrow,
        peptide_file_wide, design_file_wide, exp_file_wide, msstats_design_wide
    )
    
    for(file in files_to_check) {
        if (!file.exists(file)) stop("Required file not found: ", file)
    }
    
    return(list(
        narrow = list(
            peptide_file = peptide_file_narrow,
            design_file = design_file_narrow,
            exp_file = exp_file_narrow,
            msstats_design_file = msstats_design_narrow
        ),
        wide = list(
            peptide_file = peptide_file_wide,
            design_file = design_file_wide,
            exp_file = exp_file_wide,
            msstats_design_file = msstats_design_wide
        )
    ))
}

clean_filename <- function(x) {
  x <- gsub("\\\\", "/", x)                 # fix separators
  basename(tools::file_path_sans_ext(x))   # strip path + extension
}


#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------
# Process each experiment pair
for(i.exp in names(all.exp.pairs)) {
    message("Processing experiment pair: ", i.exp)
    
    # Get narrow and wide experiment names
    exp_narrow <- all.exp.pairs[[i.exp]]$narrow
    exp_wide <- all.exp.pairs[[i.exp]]$wide
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(exp_narrow, exp_wide)
        
        # Load narrow window data
        message("  Loading narrow window data...")
        peptide_narrow <- fread(files$narrow$peptide_file, sep = "\t", nThread = 10)
        peptide_narrow <- data.frame(peptide_narrow, check.rows = FALSE, check.names = FALSE)
        design_narrow <- read.csv(files$narrow$design_file, sep = "\t")
        exp_narrow_data <- read.csv(files$narrow$exp_file, sep = "\t")
        annot_narrow <- read.csv(files$narrow$msstats_design_file, sep = "\t")
        
        # Load wide window data  
        message("  Loading wide window data...")
        peptide_wide <- fread(files$wide$peptide_file, sep = "\t", nThread = 10)
        peptide_wide <- data.frame(peptide_wide, check.rows = FALSE, check.names = FALSE)
        design_wide <- read.csv(files$wide$design_file, sep = "\t")
        exp_wide_data <- read.csv(files$wide$exp_file, sep = "\t")
        annot_wide <- read.csv(files$wide$msstats_design_file, sep = "\t")
        
    }, error = function(e) {
        message("Error processing experiment ", i.exp, ": ", e$message)
        stop("Critical error in experiment processing - stopping execution")
    })
    
    # Process and combine data
    tryCatch({
        message("  Processing and combining datasets...")
        
        # Add batch labels to conditions and data
        annot_narrow$Condition <- paste0(annot_narrow$Condition, "_Narrow")
        annot_wide$Condition <- paste0(annot_wide$Condition, "_Wide")
        peptide_narrow$R.Condition <- paste0(peptide_narrow$R.Condition, "_Narrow")
        peptide_wide$R.Condition <- paste0(peptide_wide$R.Condition, "_Wide")
        
        # Process quantification data
        exp_narrow_data[exp_narrow_data == 0] <- NA
        exp_wide_data[exp_wide_data == 0] <- NA
        
        # Add batch suffixes to sample names
        narrow_cols <- 3:ncol(exp_narrow_data)
        wide_cols <- 3:ncol(exp_wide_data)
        colnames(exp_narrow_data)[narrow_cols] <- paste0(colnames(exp_narrow_data)[narrow_cols], "_N")
        colnames(exp_wide_data)[wide_cols] <- paste0(colnames(exp_wide_data)[wide_cols], "_W")
        
        # Log2 transformation
        exp_narrow_data[, narrow_cols] <- log2(exp_narrow_data[, narrow_cols] + 1)
        exp_wide_data[, wide_cols] <- log2(exp_wide_data[, wide_cols] + 1)
        
        # Clean metadata
        # Process narrow window data
        exp_narrow_data <- tryCatch({
            exp_narrow_data$Organism <- NULL
            row.names(exp_narrow_data) <- exp_narrow_data$Protein
            exp_narrow_data$Protein <- NULL
            exp_narrow_data
        }, error = function(e) {
            stop("Error cleaning narrow window data: ", e$message)
        })
        
        # Process wide window data
        exp_wide_data <- tryCatch({
            exp_wide_data$Organism <- NULL
            row.names(exp_wide_data) <- exp_wide_data$Protein
            exp_wide_data$Protein <- NULL
            exp_wide_data
        }, error = function(e) {
            stop("Error cleaning wide window data: ", e$message)
        })
        
        # Process design files
        design_narrow$sample_name <- paste0(design_narrow$sample_name, "_N")
        design_wide$sample_name <- paste0(design_wide$sample_name, "_W")
        design_narrow$batch <- "N"
        design_wide$batch <- "W"
        
        # Combine datasets
        exp_combined <- rbind(peptide_narrow, peptide_wide)
        annot_combined <- rbind(annot_narrow, annot_wide)
        design_combined <- rbind(design_narrow, design_wide)
        
        # Merge quantification data
        exp_mar_combined <- merge(exp_narrow_data, exp_wide_data, by = "row.names", all = TRUE)
        row.names(exp_mar_combined) <- exp_mar_combined$Row.names
        exp_mar_combined$Row.names <- NULL
        
        # Synchronize samples with design
        exp_mar_combined <- exp_mar_combined[, colnames(exp_mar_combined) %in% design_combined$sample_name]
        design_combined <- design_combined[design_combined$sample_name %in% colnames(exp_mar_combined), ]
        exp_mar_combined <- exp_mar_combined[, design_combined$sample_name]
        
    }, error = function(e) {
        message("Error in data processing: ", e$message)
        stop("Critical error in data processing - stopping execution")
    })
    
    # Generate pairwise contrasts
    conditions <- unique(design_combined$condition)
    Contrasts <- combn(conditions, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    
    message("  Processing ", length(Contrasts), " pairwise comparisons")
    
    # Process each contrast
    for(i.Contrasts in Contrasts) {
        message("    Processing contrast: ", i.Contrasts)
        
        tryCatch({
            # Parse conditions from contrast string
            i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
            i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
            
            # Count replicates per condition
            select.col.Contrasts1 <- sum(design_combined$condition == i.Contrasts1)
            select.col.Contrasts2 <- sum(design_combined$condition == i.Contrasts2)
            
            # Only process if both conditions have sufficient replicates
            if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1) {
                
                # Subset data for current contrast
                design.temp <- design_combined[design_combined$condition %in% c(i.Contrasts1, i.Contrasts2), ]
                design.temp <- design.temp[order(design.temp$condition), ]
                df.temp <- exp_mar_combined[, design.temp$sample_name]
                
                # Validate data dimensions
                if(ncol(df.temp) != nrow(design.temp)) {
                    stop("Mismatch between design matrix and expression data dimensions")
                }
                
                # Missing value imputation
                message("      Performing missing value imputation...")
                df.temp.imputed <- impute.MinDet(df.temp)
                df.temp.imputed <- data.frame(df.temp.imputed, check.rows = F, check.names = F)
                
                # Peptide count analysis (following Case study 1 approach)
                message("      Calculating peptide counts per protein...")
                pep.count.table <- tryCatch({
                    # Step 1: Count peptides per protein per file/run (as in Case study 1)
                    peptide_counts <- exp_combined %>%
                        group_by(Protein.Group, Run) %>%
                        summarise(count = n(), .groups = 'drop')
                    
                    # Handle UPS spike-ins for HEof experiments
                    if(i.exp %in% c("HEof")) {
                        peptide_counts$Protein.Group <- str_split_fixed(peptide_counts$Protein.Group, fixed("ups"), 2)[,1]
                    }
                    
                    # Filter for current contrast samples only
                    
                    peptide_counts_filtered <- peptide_counts[peptide_counts$Run %in% basename(tools::file_path_sans_ext(design.temp$file)),]
                    if(nrow(peptide_counts_filtered) == 0){
                      peptide_counts_filtered <- peptide_counts[peptide_counts$Run %in% clean_filename(design.temp$file),]
                    }
                    if(nrow(peptide_counts_filtered) == 0){
                      stop("peptide.temp is zero")
                    }
                    
             

                    # Step 2: Get minimum count per protein across all samples (as in Case study 1)
                    min_counts <- peptide_counts_filtered %>%
                        group_by(Protein.Group) %>%
                        slice_min(order_by = count, n = 1, with_ties = FALSE) %>%
                        data.frame(check.rows = FALSE, check.names = FALSE)
                    
                    # Create protein mapping for df.temp
                    df.temp.imputed$protein_ID <- str_split_fixed(row.names(df.temp.imputed), 
                                                                  fixed("|"), 3)[,2]
                    
                    
                    not.in.peptide <- df.temp.imputed[!df.temp.imputed$protein_ID %in% min_counts$Protein.Group,]
                    
                    if(nrow(not.in.peptide) != 0){
                        df.add <- data.frame(Protein.Group = not.in.peptide$protein_ID , Run = NA, count = 0)
                        min_counts <- rbind(min_counts, df.add)
                    }
                    
                    # Match proteins between expression data and peptide counts
                    min_counts_matched <- min_counts[min_counts$Protein.Group %in% df.temp.imputed$protein_ID, ]
                    df.temp.matched <- df.temp.imputed[df.temp.imputed$protein_ID %in% min_counts_matched$Protein.Group, ]
                    
                    # Ensure matching order
                    row.names(min_counts_matched) <- min_counts_matched$Protein.Group
                    min_counts_matched <- min_counts_matched[df.temp.matched$protein_ID, ]
                    
                    # Validate matching
                    if(!identical(min_counts_matched$Protein.Group, df.temp.matched$protein_ID)) {
                        stop("Protein IDs don't match between expression and peptide count data")
                    }
                    
                    # Set row names to match expression data
                    row.names(min_counts_matched) <- row.names(df.temp.matched)
                    
                    # Create final count table (add pseudocount as in Case study 1)
                    count_table <- data.frame(
                        count = min_counts_matched$count + 1,  # Add pseudocount
                        row.names = row.names(min_counts_matched)
                    )
                    
                    # Update df.temp.imputed to match filtered proteins
                    df.temp.imputed <<- df.temp.matched
                    
                    count_table
                }, error = function(e) {
                    stop("Error in peptide counting: ", e$message)
                })
                
                # Remove protein_ID column from expression data
                if("protein_ID" %in% colnames(df.temp.imputed)) {
                    df.temp.imputed$protein_ID <- NULL
                }
                
                if(!identical(row.names(df.temp.imputed) , row.names(count_table))){
                    stop("Notidentical")
                }
                
                # Statistical model setup
                message("      Setting up statistical model...")
                design_model <- 
                    model.matrix(~0 + condition + condition:batch, data = design.temp)
  
                
                colnames(design_model) <- make.names(colnames(design_model))
                
                # DEqMS Analysis
                message("      Performing DEqMS analysis...")
                results <- tryCatch({
                    # Initial limma fit
                    fit1 <- lmFit(df.temp.imputed, design = design_model)
                    
                    # Create and apply contrast
                    cont <- makeContrasts(
                        paste0("condition", i.Contrasts1, "-condition", i.Contrasts2), 
                        levels = design_model
                    )
                    fit2 <- contrasts.fit(fit1, contrasts = cont)
                    fit3 <- eBayes(fit2)
                    
                    # Add peptide count information
                    fit3$count <- count_table[rownames(fit3$coefficients), "count"]
                    
                    # Validate peptide counts
                    if(any(is.na(fit3$count)) || min(fit3$count) < 1) {
                        stop("Invalid peptide counts detected")
                    }
                    
                    # Apply DEqMS analysis
                    fit4 <- DEqMS::spectraCounteBayes(fit3)
                    
                    # Extract and format results
                    DEqMS.results <- DEqMS::outputResult(fit4, coef_col = 1)
                    
                    # Add metadata
                    DEqMS.results$contrast <- paste0(i.Contrasts1, "_vs_", i.Contrasts2)
                    DEqMS.results$analysis_date <- Sys.Date()
                    
                    DEqMS.results
                }, error = function(e) {
                    stop("Error in DEqMS analysis: ", e$message)
                })
                
                # Save results
                output_file <- file.path(OUTPUT_DIR, 
                                       paste0("DEqMS_", i.exp, "_", i.Contrasts, ".rds"))
                saveRDS(results, output_file)
                message("      Results saved to: ", output_file)
                
                # Summary statistics
                n_significant <- sum(results$adj.P.Val < 0.05, na.rm = TRUE)
                message("      Number of significant proteins (FDR < 0.05): ", n_significant)
                
                # Clean up memory
                gc()
                
            } else {
                message("      Skipping contrast due to insufficient replicates")
            }
            
        }, error = function(e) {
            message("Error in contrast analysis for ", i.Contrasts, ": ", e$message)
            stop("Critical error in DEqMS analysis - stopping execution")
        })
    }
}

message("Analysis complete!")
