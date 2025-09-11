#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Spectronaut.DEqMS.r - Differential Expression Analysis using DEqMS package
#------------------------------------------------------------------------------
# This script performs differential expression analysis on Spectronaut proteomics data
# using the DEqMS package, which accounts for peptide count information in the 
# statistical modeling of protein abundance data.
#
# Required Input Data Structure:
# - Spectronaut output directory containing:
#   * {experiment}/Spectronaut/Report.tsv (Spectronaut output)
#   * {experiment}_spt_dlfq.tsv (protein quantification data)
#   * {experiment}_spt_design.tsv (experimental design)
#   * {experiment}_spt_design_msstats.tsv (MSstats annotation)
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
    library(MSstats)
    library(DEqMS)
    library(matrixStats)
    library(multiUS)
    library(data.table)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
SPECTRONAUT_DIR <-  "../Spectronaut"  # Spectronaut directory
OUTPUT_DIR <- "Spectronaut_results/"  # Directory for results

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# List all experiment directories, excluding specified ones
all.exp <- list.files(SPECTRONAUT_DIR, pattern = "_DIA_spt_dlfq")
all.exp <- str_split_fixed(all.exp, "_DIA_spt_", 2)[,1]
all.exp <- unique(all.exp)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
    peptide_file <- file.path(SPECTRONAUT_DIR, paste0(exp_name , "_DIA"), "Spectronaut/")
    peptide_file <- file.path(peptide_file, list.files(peptide_file, pattern = "Report.tsv"))
    
    design_file <- file.path(SPECTRONAUT_DIR, paste0(exp_name, "_DIA_spt_design.tsv"))
    exp_file <- file.path(SPECTRONAUT_DIR, paste0(exp_name, "_DIA_spt_dlfq.tsv"))
    
    if (!file.exists(peptide_file)) stop("Peptide file not found: ", peptide_file)
    if (!file.exists(design_file)) stop("Design file not found: ", design_file)
    if (!file.exists(exp_file)) stop("EXp. design file not found: ", exp_file)
    
    return(list(
        peptide_file = peptide_file,
        design_file = design_file,
        exp_file = exp_file
    ))
}

clean_filename <- function(x) {
  x <- gsub("\\\\", "/", x)                 # fix separators
  basename(tools::file_path_sans_ext(x))   # strip path + extension
}

#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------

# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        peptide <- fread(files$peptide_file, sep = "\t", nThread = 15)
        peptide <- data.frame(peptide, check.rows = F, check.names = F)
        design <- read.csv(files$design_file, sep = "\t")
        exp.dlfq <- read.csv(files$exp_file, sep = "\t")
    }, error = function(e) {
        message("Error processing experiment ", i.exp, ": ", e$message)
        stop("Critical error in DEqMS analysis - stopping execution")
    })
    
    # Generate pairwise contrasts
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]

    peptide <- peptide %>% group_by(PG.ProteinAccessions, R.FileName) %>% summarise(count = n())
    peptide <- data.frame(peptide, check.rows = F, check.names = F)
    if(i.exp %in% c("HEof_w600" , "HEof_n600")){
        peptide$PG.ProteinAccessions <- str_remove(peptide$PG.ProteinAccessions, fixed("_HUMAN_UPS"))
    }
    
    # Clean sample names
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.dlfq) <- str_remove(colnames(exp.dlfq), fixed("_"))
    row.names(exp.dlfq) <- exp.dlfq$Protein
    exp.dlfq$Protein <- NULL
    exp.dlfq$Organism <- NULL
    
    exp.dlfq <- exp.dlfq[, colnames(exp.dlfq) %in% design$sample]
    design <- design[design$sample %in% colnames(exp.dlfq), ]
    
    # Process each contrast
    for(i.Contrasts in Contrasts) {
        message("  Processing contrast: ", i.Contrasts)
        groups <- str_split(i.Contrasts, "")[[1]]
        
        # Extract conditions for current contrast
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        select.col.Contrasts1 <- sum(design$condition == i.Contrasts1)
        select.col.Contrasts2 <- sum(design$condition == i.Contrasts2)
        
        
        if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1) {
        
        # Prepare design matrix for current contrast  
            # Prepare design matrix for current contrast
        design.temp <- design[design$condition %in% groups,]
        design.temp <- design.temp[order(design.temp$condition),]
        df.temp <- exp.dlfq[, design.temp$sample]
        
        peptide.temp <- peptide[peptide$R.FileName %in% basename(tools::file_path_sans_ext(design.temp$file)),]
        if(nrow(peptide.temp) == 0){
          peptide.temp <- peptide[peptide$R.FileName %in% clean_filename(design.temp$file),]
        }
        if(nrow(peptide.temp) == 0){
          stop("peptide.temp is zero")
        }
        
        peptide.temp <- peptide.temp %>% group_by(PG.ProteinAccessions) %>% 
            slice_min(order_by = count, n = 1, with_ties = FALSE)
        peptide.temp <- data.frame(peptide.temp, check.rows = F, check.names = F)
        
        # Prepare peptide count information
        
        if(i.exp %in% c( "HEof_n600" , "HEof_w600") ){
        df.temp$protein_ID <- sapply(row.names(df.temp), function(x) {
            if(grepl( "HUMAN_UPS" , toupper(x) ) | grepl( "HUMAN" , toupper(x) ) |
               grepl( "UPS" , toupper(x) )){
                x <- toupper(x)
                x <- str_split_fixed(x , fixed("|") , 3)[,3]
                x <- str_remove(x , "_HUMAN_UPS")
                x <- str_remove(x , "_HUMAN")
                x <- str_remove(x , "_UPS")
            }else{
                x <- str_split_fixed(x , fixed("|") , 3)[,2]
            } 
        }
            )
        }else{
            df.temp$protein_ID <- sapply(row.names(df.temp), function(x){
                x <- str_split(x , fixed(";"))
                x <- str_extract_all(x, "(?<=sp\\|)[^|]+")[[1]]
                x <- paste(x, collapse = ";")
            }
                
            )
        }
        
        not.in.peptide <- df.temp[!df.temp$protein_ID %in% peptide.temp$PG.ProteinAccessions,]
        
        if(nrow(not.in.peptide) != 0){
          not.in.peptide$protein_ID
          df.add <- data.frame(PG.ProteinAccessions = not.in.peptide$protein_ID , R.FileName = NA, count = 0)
          peptide.temp <- rbind(peptide.temp, df.add)
        }
        
        
        peptide.temp <- peptide.temp[peptide.temp$PG.ProteinAccessions %in% df.temp$protein_ID,]
        
        df.temp <- df.temp[df.temp$protein_ID %in% peptide.temp$PG.ProteinAccessions,]
        
        df.temp <- df.temp[!is.na(df.temp$protein_ID),]
        peptide.temp <- peptide.temp[!is.na(peptide.temp$PG.ProteinAccessions),]
        
        row.names(peptide.temp) <- peptide.temp$PG.ProteinAccessions
        peptide.temp <- peptide.temp[df.temp$protein_ID,]
        if(!identical(peptide.temp$PG.ProteinAccessions , df.temp$protein_ID)){
            stop("Not identical")
        }
        row.names(peptide.temp) <-  row.names(df.temp)
        if(!identical(row.names(peptide.temp) , row.names(df.temp))){
            stop("Not identical")
            if(!identical(peptide.temp$PG.ProteinAccessions , df.temp$protein_ID)){
                stop("Not identical")
            }
        }
        
        df.temp$protein_ID <- NULL
        
        if(nrow(peptide.temp) == 0){
            stop("peptide.temp error")
        }
        
        if(nrow(df.temp) == 0){
            stop("peptide.temp error")
        }
        
        pep.count.table <- data.frame(
            count = peptide.temp$count + 1,  # Add pseudocount
            row.names = row.names(peptide.temp)
        )
        
        df.temp[df.temp == 0] <- NA
        
        df.temp <- log2(df.temp + 1)
        
        # Impute missing values using sequential KNN
        df.temp <- impSeq(df.temp)
        df.temp <- data.frame(df.temp, check.rows = F, check.names = F)
        # Build design matrix for differential expression analysis
        design_model <- model.matrix(~0 + condition, data = design.temp)
        message("    Design matrix created with dimensions: ", 
                paste(dim(design_model), collapse = " x "))
        
        # Fit statistical models
        tryCatch({
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
            rownames(df.temp) <- df.temp$Protein.Id
            DEqMS.results$Gene.name <- df.temp[DEqMS.results$gene, ]$Gene.Names
            

            # Save results
            output_file <- file.path(OUTPUT_DIR, 
                                   paste0("DEqMS_", i.exp, "_", i.Contrasts, ".rds"))
            saveRDS(DEqMS.results, output_file)
            message("    Results saved to: ", output_file)
            
            # Clean up to free memory
            rm(fit3, fit4)
            gc()
            
        }, error = function(e) {
            message("    Error in DEqMS analysis: ", e$message)
            stop("Critical error in DEqMS analysis - stopping execution")
        })
    }
    }
    
}

#------------------------------------------------------------------------------
# End of Script
#------------------------------------------------------------------------------
message("Analysis complete!")
