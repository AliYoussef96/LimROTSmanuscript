#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# HEof.r - Multi-Method Differential Expression Analysis for Case Study 4
#------------------------------------------------------------------------------
# This script performs differential expression analysis on UPS1 Case 4 data
# using multiple statistical methods to compare different experimental conditions
# while accounting for tool and batch effects.
#
# Analysis Methods:
# - LimROTS (with trend and robust settings)
# - ROTS (standard and linear model versions)
# - Limma (with trend and robust settings)
# - SAM (Significance Analysis of Microarrays)
# - t-test
# - ANOVA
#
# Required Input:
# - UPS1.Case4 data object (loaded from package)
#
# Output:
# - RDS files containing analysis results for each method
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths and parameters
OUTPUT_DIR <- "results"
SEED <- 1597
ROTS_ITERATIONS <- 1000
SAM_PERMUTATIONS <- 100
THREADS <- 10

# Create output directory if needed
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Load required libraries
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
})

#------------------------------------------------------------------------------
# Data Loading
#------------------------------------------------------------------------------
message("Loading UPS1 Case 4 data...")
data(UPS1.Case4)

# Validate data structure
if (!all(c("tool", "Conc.", "fake.batch") %in% colnames(colData(UPS1.Case4)))) {
    stop("Missing required metadata columns in UPS1.Case4")
}

# Set random seed for reproducibility
set.seed(SEED, sample.kind = "default", kind = "default")

        
#------------------------------------------------------------------------------
# LimROTS Analysis
#------------------------------------------------------------------------------
message("Running LimROTS analysis...")


# Set metadata and formula
meta.info <- c("tool", "Conc.", "fake.batch")
K <- floor(nrow(assay(UPS1.Case4)) / 4)
group.name <- "Conc."
formula.str <- "~0+Conc.+tool+fake.batch"

# Run LimROTS analysis
UPS1_limrots <- LimROTS(
    x = UPS1.Case4,
    niter = ROTS_ITERATIONS,
    K = K,
    meta.info = meta.info,
    BPPARAM = SnowParam(THREADS, progressbar = TRUE),
    group.name = group.name,
    formula.str = formula.str,
    trend = TRUE,
    robust = TRUE,
    permutating.group = FALSE
)

saveRDS(UPS1_limrots, file.path(OUTPUT_DIR, "LimROTS.rds"))
remove(UPS1_limrots)
gc()

message("LimROTS analysis completed")


#------------------------------------------------------------------------------
# Limma Analysis
#------------------------------------------------------------------------------
message("Running limma analysis...")


# Prepare design matrix
sample_info <- colData(UPS1.Case4)
design_model <- model.matrix(~0 + Conc. + tool + fake.batch, data = sample_info)
colnames(design_model) <- make.names(colnames(design_model))

# Fit models
fit1 <- lmFit(assay(UPS1.Case4), design = design_model)
cont <- makeContrasts(contrasts = "Conc.2-Conc.1", levels = design_model)
fit2 <- contrasts.fit(fit1, contrasts = cont)
fit3 <- eBayes(fit2, trend = TRUE, robust = TRUE)

# Extract results
limma_results <- topTable(fit3, 
                        adjust = "BH", 
                        coef = "Conc.2-Conc.1",
                        sort.by = 'logFC', 
                        n = Inf)

saveRDS(limma_results, file.path(OUTPUT_DIR, "Limma.rds"))
remove(limma_results)
gc()

message("Limma analysis completed")




