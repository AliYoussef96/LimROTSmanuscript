#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Spectronaut.run.r - Multi-Method Differential Expression Analysis
#------------------------------------------------------------------------------
# This script performs differential expression analysis on Spectronaut proteomics data
# using multiple statistical methods:
# - LimROTS
# - ROTS
# - Limma
# - SAM (Significance Analysis of Microarrays)
# - t-test
# - ANOVA
#
# Required Input Data Structure:
# - Spectronaut directory containing:
#   * {experiment}_DIA_spt_dlfq.tsv (intensity data)
#   * {experiment}_DIA_spt_design.tsv (experimental design)
#
# Design file format:
# - Must contain columns: sample_name, condition, replicate
# - Conditions should be single characters (e.g., 'A', 'B', 'C')
#
# Output:
# - RDS files for each analysis method and comparison:
#   * LimROTS_{exp}_{contrast}.rds
#   * ROTS_{exp}_{contrast}.rds
#   * Limma_{exp}_{contrast}.rds
#   * SAM_{exp}_{contrast}.rds
#   * ttest_{exp}_{contrast}.rds
#   * ANOVA_{exp}_{contrast}.rds
#------------------------------------------------------------------------------

# Load required packages
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
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
INPUT_DIR <- "../Spectronaut/"  # Directory containing Spectronaut output
OUTPUT_DIR <- "Spectronaut_results/"  # Directory for results

# Analysis parameters
SEED <- 1597  # Random seed for reproducibility
ROTS_ITERATIONS <- 1000  # Number of iterations for ROTS
SAM_PERMUTATIONS <- 100  # Number of permutations for SAM
THREADS <- 10  # Number of parallel threads for LimROTS

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
    intensity_file <- file.path(INPUT_DIR, 
                               paste0(exp_name, "_DIA_spt_dlfq.tsv"))
    design_file <- file.path(INPUT_DIR, 
                            paste0(exp_name, "_DIA_spt_design.tsv"))
    
    if (!file.exists(intensity_file)) stop("Intensity file not found: ", intensity_file)
    if (!file.exists(design_file)) stop("Design file not found: ", design_file)
    
    return(list(
        intensity_file = intensity_file,
        design_file = design_file
    ))
}

#------------------------------------------------------------------------------
# Main Processing
#------------------------------------------------------------------------------
# Set random seed for reproducibility
set.seed(SEED, sample.kind = "default", kind = "default")

# List all experiments from the Spectronaut directory
all.exp <- list.files(INPUT_DIR, pattern = 'DIA_spt_dlfq')
all.exp <- str_split_fixed(all.exp, "_DIA_spt_", 2)[,1]
all.exp <- unique(all.exp)

#------------------------------------------------------------------------------
# Main Analysis
#------------------------------------------------------------------------------
# Process each experiment
for(i.exp in all.exp) {
    message("Processing experiment: ", i.exp)
    
    # Validate and read input files
    tryCatch({
        files <- validate_input_files(i.exp)
        exp.mar <- read.csv(files$intensity_file, sep = "\t")
        design <- read.csv(files$design_file, sep = "\t")
    }, error = function(e) {
        message("Error reading files for experiment ", i.exp, ": ", e$message)
        stop("Critical error in file reading - stopping execution")
    })
  
    # Preprocess expression data
    message("  Preprocessing protein expression data")
    tryCatch({
        # Set row names and clean up columns
        row.names(exp.mar) <- exp.mar$Protein
        exp.mar$Protein <- NULL
        exp.mar$Organism <- NULL
        
        # Clean sample names
        design$sample_name <- str_remove(design$sample_name, fixed("_"))
        colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
        
        # Match samples between expression and design data
        exp.mar <- exp.mar[, colnames(exp.mar) %in% design$sample_name]
        design <- design[design$sample_name %in% colnames(exp.mar),]
        exp.mar <- exp.mar[, design$sample_name]
    }, error = function(e) {
        message("  Error in data preprocessing: ", e$message)
        stop("Critical error in data preprocessing - stopping execution")
    })

    # Generate pairwise contrasts
    message("  Generating condition contrasts")
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    message("  Processing ", length(Contrasts), " pairwise comparisons")
    
    # Process each contrast
    for(i.Contrasts in Contrasts) {
        message("  Processing contrast: ", i.Contrasts)
        
        # Extract conditions for current contrast
        groups <- str_split(i.Contrasts, "")[[1]]
        i.Contrasts1 <- groups[1]
        i.Contrasts2 <- groups[2]
        
        # Count samples per group
        select.col.Contrasts1 <- sum(design$condition == i.Contrasts1)
        select.col.Contrasts2 <- sum(design$condition == i.Contrasts2)
        
        # Only process if both conditions have sufficient replicates
        if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1) {
            message("    Found ", select.col.Contrasts1, " samples in group ", 
                   i.Contrasts1, " and ", select.col.Contrasts2, 
                   " samples in group ", i.Contrasts2)
            # Subset data for current contrast
            design.temp <- design[design$condition %in% groups,]
            df.temp <- exp.mar[, design.temp$sample_name]
            
            # Data preprocessing
            message("    Preprocessing data")
            tryCatch({
                # Convert zeros to NA and log2 transform
                df.temp[df.temp == 0] <- NA
                df.temp <- log2(df.temp + 1)
                # Impute missing values
                df.temp <- impSeq(df.temp)
                df.temp <- data.frame(df.temp, check.names = FALSE)
                
                # Create sample metadata
                sample_info <- data.frame(
                    sample_id = colnames(df.temp),
                    group = factor(c(rep(i.Contrasts1, select.col.Contrasts1),
                                   rep(i.Contrasts2, select.col.Contrasts2)))
                )
                rownames(sample_info) <- sample_info$sample_id
                
                # Create feature metadata
                feature_info <- data.frame(
                    protein_id = rownames(df.temp)
                )
                rownames(feature_info) <- feature_info$protein_id
            }, error = function(e) {
                message("    Error in data preprocessing: ", e$message)
                return(NULL)
            })
      
            # Create SummarizedExperiment object
            message("    Running statistical analyses")
            
            # 1. LimROTS Analysis
            message("      Running LimROTS analysis")
            tryCatch({
                se <- SummarizedExperiment(
                    assays = list(protein.exp = df.temp),
                    colData = sample_info,
                    rowData = feature_info 
                )
                
                se <- LimROTS(
                    x = se,
                    niter = ROTS_ITERATIONS,
                    K = floor(nrow(df.temp) / 2),
                    meta.info = "group",
                    BPPARAM = MulticoreParam(THREADS, progressbar = TRUE),
                    group.name = "group",
                    formula.str = "~ 0 + group",
                    trend = TRUE,
                    robust = TRUE
                )
                saveRDS(se, file.path(OUTPUT_DIR, sprintf("LimROTS_%s_%s.rds", i.exp, i.Contrasts)))
                rm(se); gc()
            }, error = function(e) {
                message("      Error in LimROTS analysis: ", e$message)
            })
            
            # 2. ROTS Analysis
            message("      Running ROTS analysis")
            tryCatch({
                rots_results <- ROTS(
                    data = df.temp,
                    groups = as.numeric(sample_info$group),
                    B = ROTS_ITERATIONS,
                    K = floor(nrow(df.temp) / 2),
                    seed = SEED,
                    progress = TRUE,
                    verbose = TRUE
                )
                saveRDS(rots_results, file.path(OUTPUT_DIR, sprintf("ROTS_%s_%s.rds", i.exp, i.Contrasts)))
                rm(rots_results); gc()
            }, error = function(e) {
                message("      Error in ROTS analysis: ", e$message)
            })
            
            # 3. Limma Analysis
            message("      Running Limma analysis")
            tryCatch({
                design_matrix <- model.matrix(~0 + group, data = sample_info)
                colnames(design_matrix) <- levels(sample_info$group)
                fit <- lmFit(df.temp, design = design_matrix)
                contrast_matrix <- makeContrasts(
                    contrasts = paste0(i.Contrasts2, "-", i.Contrasts1),
                    levels = design_matrix
                )
                fit <- eBayes(contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)
                limma_results <- topTable(fit, adjust = "BH", sort.by = "logFC", n = Inf)
                saveRDS(limma_results, file.path(OUTPUT_DIR, sprintf("Limma_%s_%s.rds", i.exp, i.Contrasts)))
                rm(limma_results, fit); gc()
            }, error = function(e) {
                message("      Error in Limma analysis: ", e$message)
            })
            
            # 4. SAM Analysis
            message("      Running SAM analysis")
            tryCatch({
                sam_data <- list(
                    x = as.matrix(df.temp),
                    y = as.numeric(sample_info$group),
                    geneid = feature_info$protein_id,
                    genenames = feature_info$protein_id,
                    logged2 = TRUE
                )
                sam_obj <- samr(sam_data, resp.type = "Two class unpaired", nperms = SAM_PERMUTATIONS)
                logFC <- log2(sam_obj$foldchange)
                pval <- samr.pvalues.from.perms(sam_obj$tt, sam_obj$ttstar)
                adj_p <- p.adjust(pval, method = "BH")
                sam_results <- data.frame(
                    logFC = logFC,
                    pvalue = pval,
                    adj.pvalue = adj_p,
                    row.names = feature_info$protein_id
                )
                saveRDS(sam_results, file.path(OUTPUT_DIR, sprintf("SAM_%s_%s.rds", i.exp, i.Contrasts)))
                rm(sam_results, sam_obj, sam_data); gc()
            }, error = function(e) {
                message("      Error in SAM analysis: ", e$message)
            })
            
            # 5. t-test Analysis
            message("      Running t-test analysis")
            tryCatch({
                ttest_results <- t(sapply(1:nrow(df.temp), function(i) {
                    test <- t.test(as.numeric(df.temp[i,]) ~ sample_info$group)
                    fc <- mean(df.temp[i, sample_info$group == i.Contrasts2]) - 
                         mean(df.temp[i, sample_info$group == i.Contrasts1])
                    c(pval = test$p.value, fc = fc)
                }))
                ttest_results <- as.data.frame(ttest_results)
                ttest_results$adj.pvalue <- p.adjust(ttest_results$pval, method = "BH")
                rownames(ttest_results) <- rownames(df.temp)
                saveRDS(ttest_results, file.path(OUTPUT_DIR, sprintf("ttest_%s_%s.rds", i.exp, i.Contrasts)))
                rm(ttest_results); gc()
            }, error = function(e) {
                message("      Error in t-test analysis: ", e$message)
            })
            
            # 6. ANOVA Analysis
            message("      Running ANOVA analysis")
            tryCatch({
                anova_results <- do.call(rbind, lapply(1:nrow(df.temp), function(i) {
                    test_data <- df.temp[i,]
                    fit <- aov(as.numeric(test_data) ~ group, data = sample_info)
                    summ <- summary(fit)[[1]]
                    fc <- mean(test_data[sample_info$group == i.Contrasts2]) - 
                         mean(test_data[sample_info$group == i.Contrasts1])
                    data.frame(
                        logFC = fc,
                        P.Value = summ$`Pr(>F)`[1],
                        row.names = rownames(df.temp)[i]
                    )
                }))
                anova_results$adj.P.Val <- p.adjust(anova_results$P.Value, method = "BH")
                saveRDS(anova_results, file.path(OUTPUT_DIR, sprintf("ANOVA_%s_%s.rds", i.exp, i.Contrasts)))
                rm(anova_results); gc()
            }, error = function(e) {
                message("      Error in ANOVA analysis: ", e$message)
            })
            
            message("    Completed all analyses")
        } else {
            message("    Skipping contrast - insufficient replicates")
        }
        message("  Completed contrast: ", i.Contrasts)
    }
    message("Completed experiment: ", i.exp)
}
message("All analyses completed successfully!")
