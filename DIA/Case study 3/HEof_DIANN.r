#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# HEof_DIANN.r - Multi-Method Differential Expression Analysis for HEof Study
#------------------------------------------------------------------------------
# This script performs differential expression analysis comparing narrow and wide
# window DIA-NN proteomics data using multiple statistical methods:
# - LimROTS
# - ROTS
# - Limma
# - SAM (Significance Analysis of Microarrays)
#
# Required Input Data Structure:
# - HEof_n600_DIA_DIANN_dlfq.tsv (Narrow window data)
# - HEof_w600_DIA_DIANN_dlfq.tsv (Wide window data)
# - HEof_*_DIA_DIANN_design.tsv (Experimental designs)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Window-specific suffixes will be added (_N for narrow, _W for wide)
#
# Output:
# - RDS files for each analysis method and comparison
# - MDS plots for quality control
#------------------------------------------------------------------------------

# Load required packages
suppressPackageStartupMessages({
    library(limma)
    library(stringr)
    library(SummarizedExperiment)
    library(LimROTS)
    library(ROTS)
    library(BiocParallel)
    library(rrcovNA)
    library(imputeLCMD)
    library(samr)
    library(DEP)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths and parameters
DIANN_DIR <-  "../DIANN"
OUTPUT_DIR <- "DIANN_results"

# Analysis parameters
SEED <- 1597
ROTS_ITERATIONS <- 1000
SAM_PERMUTATIONS <- 100
THREADS <- 10

# Create output directories
for (dir in c(OUTPUT_DIR)) {
    if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}

# Set random seed for reproducibility
set.seed(SEED, sample.kind = "default", kind = "default")

#------------------------------------------------------------------------------
# Data Loading
#------------------------------------------------------------------------------
message("Loading DIA-NN output data...")

# Load narrow window data
exp.mar1 <- tryCatch({
    data <- read.csv(file.path(DIANN_DIR, "HEof_n600_DIA_DIANN_dlfq.tsv"), sep = "\t")
    data[data == 0] <- NA  # Replace 0s with NA
    data
}, error = function(e) {
    stop("Error loading narrow window data: ", e$message)
})

# Load wide window data
exp.mar2 <- tryCatch({
    data <- read.csv(file.path(DIANN_DIR, "HEof_w600_DIA_DIANN_dlfq.tsv"), sep = "\t")
    data[data == 0] <- NA  # Replace 0s with NA
    data
}, error = function(e) {
    stop("Error loading wide window data: ", e$message)
})

#------------------------------------------------------------------------------
# Data Preprocessing
#------------------------------------------------------------------------------
message("Preprocessing data...")

# Add window-specific suffixes
message("Adding window identifiers...")
colnames(exp.mar1)[3:22] <- paste0(colnames(exp.mar1)[3:22], "_N")
colnames(exp.mar2)[3:26] <- paste0(colnames(exp.mar2)[3:26], "_W")

# Log2 transformation
message("Performing log2 transformation...")
exp.mar1[, 3:22] <- log2(exp.mar1[, 3:22] + 1)
exp.mar2[, 3:26] <- log2(exp.mar2[, 3:26] + 1)


#------------------------------------------------------------------------------
# Data Cleanup
#------------------------------------------------------------------------------
message("Cleaning and formatting data...")

# Process narrow window data
exp.mar1 <- tryCatch({
    exp.mar1$Organism <- NULL
    row.names(exp.mar1) <- exp.mar1$Protein
    exp.mar1$Protein <- NULL
    exp.mar1
}, error = function(e) {
    stop("Error cleaning narrow window data: ", e$message)
})

# Process wide window data
exp.mar2 <- tryCatch({
    exp.mar2$Organism <- NULL
    row.names(exp.mar2) <- exp.mar2$Protein
    exp.mar2$Protein <- NULL
    exp.mar2
}, error = function(e) {
    stop("Error cleaning wide window data: ", e$message)
})

#------------------------------------------------------------------------------
# Design File Processing
#------------------------------------------------------------------------------
message("Processing experimental design files...")

# Load design files
design1 <- tryCatch({
    design <- read.csv(file.path(DIANN_DIR, "HEof_n600_DIA_DIANN_design.tsv"), sep = "\t")
    if (!all(c("sample_name", "condition", "replicate") %in% colnames(design))) {
        stop("Missing required columns in narrow window design file")
    }
    design
}, error = function(e) {
    stop("Error loading narrow window design: ", e$message)
})

design2 <- tryCatch({
    design <- read.csv(file.path(DIANN_DIR, "HEof_w600_DIA_DIANN_design.tsv"), sep = "\t")
    if (!all(c("sample_name", "condition", "replicate") %in% colnames(design))) {
        stop("Missing required columns in wide window design file")
    }
    design
}, error = function(e) {
    stop("Error loading wide window design: ", e$message)
})

# Add window identifiers
design1$sample_name <- paste0(design1$sample_name, "_N")
design2$sample_name <- paste0(design2$sample_name, "_W")
design1$batch <- "N"
design2$batch <- "W"

# Combine designs
design <- rbind(design1, design2)

#------------------------------------------------------------------------------
# Data Integration
#------------------------------------------------------------------------------
message("Integrating datasets...")

# Merge expression data
exp.mar <- tryCatch({
    # Merge matrices
    merged <- merge(exp.mar1, exp.mar2, by = "row.names", all = TRUE)
    row.names(merged) <- merged$Row.names
    merged$Row.names <- NULL
    
    # Synchronize with design
    merged <- merged[, colnames(merged) %in% design$sample_name]
    design <- design[design$sample_name %in% colnames(merged), ]
    merged <- merged[, design$sample_name]
    
    merged
}, error = function(e) {
    stop("Error merging datasets: ", e$message)
})


###########

#------------------------------------------------------------------------------
# Setup for Differential Expression Analysis
#------------------------------------------------------------------------------
message("Setting up differential expression analysis...")

# Generate pairwise contrasts
Contrasts <- combn(unique(design$condition), 2, simplify = FALSE)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
Contrasts <- Contrasts[!sapply(Contrasts, function(x) {
    length(unique(strsplit(x, "")[[1]])) == 1
})]

# Final data synchronization
exp.mar <- exp.mar[, colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar), ]
exp.mar <- exp.mar[, design$sample_name]

i.exp <- "HEof"

#------------------------------------------------------------------------------
# Differential Expression Analysis
#------------------------------------------------------------------------------
message("Beginning differential expression analysis...")

for(i.Contrasts in Contrasts) {
    message(sprintf("Processing contrast: %s", i.Contrasts))
    
        # Parse contrast conditions
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        # Check sample sizes
        select.col.Contrasts1 <- nrow(design[design$condition %in% c(i.Contrasts1), ])
        select.col.Contrasts2 <- nrow(design[design$condition %in% c(i.Contrasts2), ])
        
        if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1) {
            message("Sufficient samples found for both conditions")
            
            # Prepare data for current contrast
            design.temp <- design[design$condition %in% c(i.Contrasts1, i.Contrasts2), ]
            df.temp <- exp.mar[, colnames(exp.mar) %in% design.temp$sample_name]
            df.temp <- df.temp[, design.temp$sample_name]
            
            # Impute missing values
            message("Imputing missing values...")
            df.temp <- impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = FALSE, check.rows = FALSE)

            # Add random effect to a subset of ECOLI proteins in 'W' batch
            add.effect <- design.temp[design.temp$batch == "W" & design.temp$condition == design.temp$condition[1],]
            ecoli <- sample(which(grepl("ECOLI", row.names(df.temp))), 500)
            df.temp[ecoli,which(colnames(df.temp) %in% add.effect$sample_name)] <- 
                df.temp[ecoli,which(colnames(df.temp) %in% add.effect$sample_name)] + runif(500, min = 5, max = 20)
            row.names(df.temp)[ecoli] <- paste0(row.names(df.temp)[ecoli], "_added")
        
        # Prepare sample metadata
        sample_info <- data.frame(
            sample_id = colnames(df.temp),
            group = design.temp$condition,
            batches = design.temp$batch)
        
        rownames(sample_info) <- row.names(sample_info$sample_id)
        sample_info$group <- as.factor( as.numeric(as.factor(sample_info$group)) )
        sample_info$batches <- as.factor(sample_info$batches)
        
        # Prepare feature metadata
        feature_info <- data.frame(
            protein_id = row.names(df.temp)
        )
        rownames(feature_info) <- feature_info$protein_id
        
        # run ANOVA
        anova.results <- data.frame()
        for(i.test in seq_len(nrow(df.temp))){
          i.test.df <- df.temp[i.test,]
          fit = aov(as.numeric(i.test.df)~group+group:batches, data = sample_info)
          fit <- summary(fit)[[1]]
          fc.calc <- mean( as.numeric( i.test.df[,1:select.col.Contrasts1] ) ) - mean( as.numeric(i.test.df[,select.col.Contrasts1+1:select.col.Contrasts2]) )
          anova.results = rbind(anova.results, data.frame(row.names = row.names(i.test.df),
                                                          logFC=fc.calc,
                                                          P.Value=fit$`Pr(>F)`[1]) ) }
        anova.results$adj.P.Val <- p.adjust(anova.results$P.Value, method  = "BH")
        saveRDS(anova.results, paste0("DIANN_results/", "ANOVA_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(anova.results)
        gc()
        # run DEP



 
    } else {
        warning(sprintf("Insufficient samples for contrast %s: %d vs %d samples", 
                       i.Contrasts, select.col.Contrasts1, select.col.Contrasts2))
    }

}

message("Differential expression analysis completed.")
