
#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# FragPipe.run.r - Multi-Method Differential Expression Analysis
#------------------------------------------------------------------------------
# This script performs differential expression analysis on FragPipe proteomics data
# using multiple statistical methods:
# - LimROTS
# - ROTS
# - Limma
# - SAM (Significance Analysis of Microarrays)
# - t-test
# - ANOVA
#
# Required Input Data Structure:
# - FragPipe directory containing:
#   * {experiment}_LFQ_FragPipe_dlfq_pro_intensity.tsv (intensity data)
#   * {experiment}_FragPipe_design.tsv (experimental design)
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
    library(multiUS)
})

#------------------------------------------------------------------------------
# Configuration
#------------------------------------------------------------------------------
# Set paths (modify these as needed)
INPUT_DIR <- "FragPipe/"  # Directory containing FragPipe output
OUTPUT_DIR <- "FragPipe_results/"  # Directory for results

# Analysis parameters
SEED <- 1597  # Random seed for reproducibility
ROTS_ITERATIONS <- 1000  # Number of iterations for ROTS
SAM_PERMUTATIONS <- 100  # Number of permutations for SAM
THREADS <- 10  # Number of parallel threads for LimROTS

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Helper Functions
#------------------------------------------------------------------------------
# Function to validate input files
validate_input_files <- function(exp_name) {
    intensity_file <- file.path(INPUT_DIR, 
                               paste0(exp_name, "_LFQ_FragPipe_dlfq_pro_intensity.tsv"))
    design_file <- file.path(INPUT_DIR, 
                            paste0(exp_name, "_LFQ_FragPipe_design.tsv"))
    
    if (!file.exists(intensity_file)) stop("Intensity file not found: ", intensity_file)
    if (!file.exists(design_file)) stop("Design file not found: ", design_file)
    
    return(list(
        intensity_file = intensity_file,
        design_file = design_file
    ))
}

# Set random seed for reproducibility
set.seed(SEED, sample.kind = "default", kind = "default")

# List all experiments from the FragPipe directory
all.exp <- list.files(INPUT_DIR,  pattern = '_LFQ_FragPipe_dlfq_pro_intensity')
all.exp <- str_split_fixed(all.exp, "_LFQ_FragPipe_", 2)[,1]
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
  
    # Generate pairwise contrasts
    message("  Generating condition contrasts")
    Contrasts <- combn(design$condition, 2, simplify = FALSE)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    # Remove self-comparisons
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) 
        length(unique(strsplit(x, "")[[1]])) == 1)]
    message("  Processing ", length(Contrasts), " pairwise comparisons")
    
    # Preprocess expression data
    message("  Preprocessing protein expression data")
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    
    # Clean sample names
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Ensure data and design matrix match
    exp.mar <- exp.mar[, colnames(exp.mar) %in% design$sample_name]
    design <- design[design$sample_name %in% colnames(exp.mar),]
    exp.mar <- exp.mar[, design$sample_name]
  
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
                df.temp <- seqKNNimp(df.temp)
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
                stop("Critical error in data preprocessing - stopping execution")
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
                saveRDS(se, file.path(OUTPUT_DIR, 
                                    paste0("LimROTS_", i.exp, "_", i.Contrasts, ".rds")))
                rm(se); gc()
            }, error = function(e) {
                message("      Error in LimROTS analysis: ", e$message)
                stop("Critical error in LimROTS analysis - stopping execution")
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
                saveRDS(rots_results, file.path(OUTPUT_DIR, 
                                              paste0("ROTS_", i.exp, "_", i.Contrasts, ".rds")))
                rm(rots_results); gc()
            }, error = function(e) {
                message("      Error in ROTS analysis: ", e$message)
                stop("Critical error in ROTS analysis - stopping execution")
            })
            
            # 3. Limma Analysis
            message("      Running Limma analysis")
            tryCatch({
                design.model <- model.matrix(~0 + group, data = sample_info)
                fit <- lmFit(df.temp, design = design.model)
                contrast <- makeContrasts(
                    contrasts = paste0("group", i.Contrasts1, "-group", i.Contrasts2),
                    levels = design.model
                )
                fit <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
                limma.res <- topTable(fit, adjust = "BH", sort.by = "logFC", n = Inf)
                saveRDS(limma.res, file.path(OUTPUT_DIR, 
                                           paste0("Limma_", i.exp, "_", i.Contrasts, ".rds")))
                rm(limma.res); gc()
            }, error = function(e) {
                message("      Error in Limma analysis: ", e$message)
                stop("Critical error in Limma analysis - stopping execution")
            })
            
            # 4. SAM Analysis
            message("      Running SAM analysis")
            tryCatch({
                sam.data <- list(
                    x = as.matrix(df.temp),
                    y = as.numeric(sample_info$group),
                    geneid = feature_info$protein_id,
                    genenames = feature_info$protein_id,
                    logged2 = TRUE
                )
                sam.obj <- samr(sam.data, resp.type = "Two class unpaired", 
                               nperms = SAM_PERMUTATIONS)
                logFC <- log2(sam.obj$foldchange)
                pval <- samr.pvalues.from.perms(sam.obj$tt, sam.obj$ttstar)
                adj.p <- p.adjust(pval, method = "BH")
                SAM.res <- cbind(logFC, pvalue = pval, adj.pvalue = adj.p)
                saveRDS(SAM.res, file.path(OUTPUT_DIR, 
                                         paste0("SAM_", i.exp, "_", i.Contrasts, ".rds")))
                rm(SAM.res); gc()
            }, error = function(e) {
                message("      Error in SAM analysis: ", e$message)
            })
            
            # 5. t-test Analysis
            message("      Running t-test analysis")
            tryCatch({
                ttest.res <- t(sapply(seq_len(nrow(df.temp)), function(i){
                    test <- t.test(as.numeric(df.temp[i,]) ~ sample_info$group)
                    fc <- mean(df.temp[i, sample_info$group == i.Contrasts1]) - 
                          mean(df.temp[i, sample_info$group == i.Contrasts2])
                    c(pval = test$p.value, fc = fc)
                }))
                ttest.res <- as.data.frame(ttest.res)
                ttest.res$adj.pvalue <- p.adjust(ttest.res$pval, method = "BH")
                rownames(ttest.res) <- rownames(df.temp)
                saveRDS(ttest.res, file.path(OUTPUT_DIR, 
                                           paste0("ttest_", i.exp, "_", i.Contrasts, ".rds")))
                rm(ttest.res); gc()
            }, error = function(e) {
                message("      Error in t-test analysis: ", e$message)
                stop("Critical error in t-test analysis - stopping execution")
            })
            
            # 6. ANOVA Analysis
            message("      Running ANOVA analysis")
            tryCatch({
                anova.results <- data.frame()
                for(i.test in seq_len(nrow(df.temp))) {
                    i.test.df <- df.temp[i.test,]
                    fit <- aov(as.numeric(i.test.df) ~ group, data = sample_info)
                    fit <- summary(fit)[[1]]
                    fc.calc <- mean(as.numeric(i.test.df[,1:select.col.Contrasts1])) - 
                              mean(as.numeric(i.test.df[,(select.col.Contrasts1+1):
                                                      (select.col.Contrasts1+select.col.Contrasts2)]))
                    anova.results <- rbind(anova.results, 
                                         data.frame(row.names = row.names(i.test.df),
                                                  logFC = fc.calc,
                                                  P.Value = fit$`Pr(>F)`[1]))
                }
                anova.results$adj.P.Val <- p.adjust(anova.results$P.Value, method = "BH")
                saveRDS(anova.results, file.path(OUTPUT_DIR, 
                                               paste0("ANOVA_", i.exp, "_", i.Contrasts, ".rds")))
                rm(anova.results); gc()
            }, error = function(e) {
                message("      Error in ANOVA analysis: ", e$message)
                stop("Critical error in ANOVA analysis - stopping execution")
            })
            
            message("    Completed all analyses for contrast ", i.Contrasts)
        } else {
            message("    Skipping contrast due to insufficient replicates")
        }
    }
}

message("Analysis complete!")
