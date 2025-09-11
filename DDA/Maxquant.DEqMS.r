#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Maxquant.DEqMS.r - Differential Expression Analysis using DEqMS package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on Maxquant proteomics data
# using the DEqMS package, which accounts for peptide count information in the 
# statistical modeling of protein abundance data.
#
# Required Input Data Structure:
# - Maxquant output directory containing:
#   * {experiment}/Maxquant/TOP0/combined_protein.tsv (protein quantification data)
#   * {experiment}_Maxquant_design.tsv (experimental design)
#   * {experiment}_Maxquant_design_msstats.tsv (MSstats annotation)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions should be single characters (e.g., 'A', 'B', 'C')
#
# Output:
# - RDS files containing DEqMS analysis results for each comparison
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
    library(DEqMS)
    library(matrixStats)
    library(multiUS)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
INPUT_DIR <- "Maxquant"  # Base directory for experiments
OUTPUT_DIR <- "Maxquant_results/"  # Directory for results

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
    peptide_file <- file.path(INPUT_DIR,  paste0(exp_name , "_LFQ") , "Maxquant/proteinGroups.txt")
    design_file <- file.path(INPUT_DIR, paste0(exp_name, "_LFQ_Maxquant_design.tsv"))
    exp_file <- file.path(INPUT_DIR, paste0(exp_name, "_LFQ_Maxquant_dlfq_pro_intensity.tsv"))
    
    if (!file.exists(peptide_file)) stop("Protein file not found: ", peptide_file)
    if (!file.exists(design_file)) stop("Design file not found: ", design_file)
    if (!file.exists(exp_file)) stop("Exp. file not found: ", design_file)
    
    return(list(
        peptide_file = peptide_file,
        design_file = design_file,
        exp_file = exp_file
    ))
}

# List all experiment folders except specific ones
all.exp <- list.files(INPUT_DIR, pattern = '_LFQ_Maxquant_dlfq_pro_intensity')
all.exp <- unique(str_split_fixed(all.exp, "_LFQ_Maxquan", 2)[,1])

#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------
# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        peptide <- read.csv(files$peptide_file, sep = "\t")
        design <- read.csv(files$design_file, sep = "\t")
        exp.dlfq <- read.csv(files$exp_file, sep = "\t")
    }, error = function(e) {
        message("Error processing experiment ", i.exp, ": ", e$message)
        stop("Critical error in experiment processing - stopping execution")
    })
    
    # Generate pairwise contrasts
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    
    # Format protein data
    peptide <- peptide[!peptide$Reverse == "+",]
    peptide <- peptide[!peptide$Potential.contaminant == "+",]
    
    row.names(peptide) <- peptide$Protein
    peptide$Protein <- NULL
    peptide$Organism <- NULL
    
    row.names(exp.dlfq) <- exp.dlfq$Protein
    exp.dlfq$Protein <- NULL
    exp.dlfq$Organism <- NULL
    
    # Clean sample names
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(peptide) <- str_remove(colnames(peptide), fixed("_"))
    
    # Process each contrast
    for(i.Contrasts in Contrasts) {
        message("  Processing contrast: ", i.Contrasts)
        
        groups <- str_split(i.Contrasts, "")[[1]]
        # Extract conditions for current contrast
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        select.col.Contrasts1 <- sum(design$condition == i.Contrasts1)
        select.col.Contrasts2 <- sum(design$condition == i.Contrasts2)
        
        # Prepare design matrix for current contrast
        if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1) {
            
            tryCatch({
                design.temp <- design[design$condition %in% groups,]
                df.temp <- exp.dlfq[, design.temp$sample_name]
                
                message("Number of samples : ", ncol(df.temp))
                
                if(ncol(df.temp) <= 3){
                    stop("df.temp Less than 3!!")
                }
                
                peptide.razor <- which( grepl("Razor...unique.peptides.",colnames(peptide)) )
                colnames(peptide)[peptide.razor]
                
                matches <- str_remove(colnames(peptide)[peptide.razor] , fixed("Razor...unique.peptides.") )
                matches <- matches[matches %in% design.temp$sample_name]
                if(length(matches) == 0){
                  stop("Match is zero")
                }
                matches <- paste0("Razor...unique.peptides." , matches )
   
                # Prepare peptide count information
                pep.count.table <- data.frame(
                  count = rowMins(as.matrix(peptide[,matches])) + 1,  # Add pseudocount
                  row.names = row.names(peptide),
                  proteins = row.names(peptide)
                )
                
                not.in.peptide <- df.temp[!row.names(df.temp) %in% pep.count.table$proteins,]
                
                if(nrow(not.in.peptide) != 0){
                    df.add <- data.frame(proteins = row.names(not.in.peptide) , count = 0)
                    pep.count.table <- rbind(pep.count.table, df.add)
                }
                
                pep.count.table <- pep.count.table[row.names(pep.count.table) %in% row.names(df.temp),]
                pep.count.table <- pep.count.table[!is.na(pep.count.table$proteins),]
                
                
                df.temp <- df.temp[row.names(df.temp) %in% row.names(pep.count.table),]
                df.temp <- na.omit(df.temp)
                if(nrow(df.temp) == 0){
                  stop("df.temp is Zero")
                }
                
                pep.count.table <- pep.count.table[row.names(df.temp),]
                pep.count.table$proteins <- NULL
                
                if(!identical(row.names(pep.count.table) , row.names(df.temp))){
                    stop("Not Identical")
                }
                
                # Data preprocessing
                # Log2 transform LFQ intensities with pseudocount
                df.temp[df.temp == 0] <- NA
                df.temp <- log2(df.temp + 1)
                # Impute missing values
                df.temp <- impSeq(df.temp)
                df.temp <- data.frame(df.temp, check.names = FALSE)
                
                # Build design matrix for differential expression analysis
                design_model <- model.matrix(~0 + condition, data = design.temp)
                message("    Design matrix created with dimensions: ", 
                        paste(dim(design_model), collapse = " x "))
                
                # Fit statistical models
                
                # Initial limma fit
                fit1 <- lmFit(df.temp, design = design_model)
                
                # Create and apply contrast
                cont <- makeContrasts(
                    paste0("condition", i.Contrasts1, "-condition", i.Contrasts2), 
                    levels = design_model
                )
                fit2 <- contrasts.fit(fit1, contrasts = cont)
                fit3 <- eBayes(fit2)
                
                # Add peptide count information
                fit3$count <- pep.count.table[rownames(fit3$coefficients), "count"]
                
                # Validate peptide counts
                if(min(fit3$count) < 1) {
                    stop("Invalid peptide counts detected (< 1)")
                }
                
                # Apply DEqMS analysis
                fit4 <- spectraCounteBayes(fit3)
                
                # Extract and format results
                DEqMS.results <- outputResult(fit4, coef_col = 1)
                rownames(df.temp) <- df.temp$Majority.protein.IDs
                DEqMS.results$Gene.name <- df.temp[DEqMS.results$gene, ]$Gene.names
                
                # Save results
                output_file <- file.path(OUTPUT_DIR, 
                                         paste0("DEqMS_", i.exp, "_", i.Contrasts, ".rds"))
                saveRDS(DEqMS.results, output_file)
                message("    Results saved to: ", output_file)
                
                # Clean up to free memory
                rm(fit3, fit4)
                gc()
            }, error = function(e) {
                message("Error in DEqMS analysis for contrast ", i.Contrasts, ": ", e$message)
                stop("Critical error in DEqMS analysis - stopping execution")
            })
            
        } else {
            message("    Skipping contrast due to missing samples in LFQ data")
        }
    }
}



message("Analysis complete!")
