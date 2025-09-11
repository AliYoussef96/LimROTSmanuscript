# Configuration
#------------------------------------------------------------------------------

# Set paths based on analysis mode
RESULTS_DIR <- "results"
OUTPUT_DIR <- "metrics"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Load required libraries
suppressPackageStartupMessages({
    library(stringr)  
    library(pROC)    
    library(caret) 
})

# Ground truth function for DIA datasets
get_ground_truth_DIA <- function(dataset_name) {
  # Ground truth mapping for DIA datasets based on trueDE column
  ground_truth_map <- list(
    # Human+Yeast+E.coli datasets - trueDE: YEAST;ECOLI
    "HYEtims735" = c("YEAST", "ECOLI"),
    
    # Mouse+Yeast datasets - trueDE: MOUSE
    "MYtims709" = c("MOUSE"),
    
    # Human+E.coli+UPS datasets - trueDE: UPS (both narrow and wide)
    "HEof" = c("UPS"),

    # Human+Yeast datasets - trueDE: YEAST
    "HYtims134" = c("YEAST"),
    
    # Human+E.coli datasets - trueDE: ECOLI
    "HEqe777" = c("ECOLI"),
    "HEqe408" = c("ECOLI")
  )
  
  if (dataset_name %in% names(ground_truth_map)) {
    return(ground_truth_map[[dataset_name]])
  } else {
    warning(paste("Dataset", dataset_name, "not found in ground truth mapping"))
    return(NULL)
  }
}



#------------------------------------------------------------------------------
# List differential expression result files
results.files <- list.files(RESULTS_DIR)

# Set FDR threshold for classification
fdr_threshold <- 0.05

# Initialize dataframe to store performance metrics
metrics.all <- data.frame(row.names = 1:17)

# Loop over each result file
for(i.file in results.files){
    # Extract tool, experiment, and contrast from filename
    


    dea.tool <- str_remove(i.file , ".rds")
    
    # Get ground truth organisms for this experiment
    ground_truth_organisms <- 'UPS'
    
    if(is.null(ground_truth_organisms)) {
        warning(paste("Skipping", i.file, "- no ground truth found for", exp.name))
        stop("Critical error in ground truth mapping - stopping execution")
    }
    
    # Load DEA results
    result.read <- readRDS(file.path(RESULTS_DIR, i.file))
    
    if(dea.tool == "MSstats"){
    result.read$Protein <- str_split_fixed(result.read$Protein, fixed(";"), 2)[,1]
    result.read.see <- result.read[!is.na(result.read$pvalue),]
    if(nrow(result.read.see) != 0){
        result.read$issue <- ifelse(is.na(result.read$issue) , "OK" , result.read$issue ) 
        result.read <- result.read[result.read$issue != "completeMissing",]  
        result.read$pvalue <- ifelse(is.na(result.read$pvalue) & result.read$issue == "oneConditionMissing"
                                     ,min(result.read$pvalue, na.rm = T) , result.read$pvalue )
        result.read$pvalue <- ifelse(is.na(result.read$pvalue) 
                                     ,max(result.read$pvalue, na.rm = T) , result.read$pvalue )
    }
    }
    
    # Determine ground truth based on protein identifiers
    if(length(ground_truth_organisms) == 2){
        # Two organisms expected to be differentially expressed
        if(dea.tool == "ROTS"){
            ground_truth <- ifelse(grepl(ground_truth_organisms[1], rownames(result.read$data), ignore.case = TRUE) | 
                                 grepl(ground_truth_organisms[2], rownames(result.read$data), ignore.case = TRUE), 1, 0)
        } else if(dea.tool == "MSstats"){
            ground_truth <- ifelse(grepl(ground_truth_organisms[1], result.read$Protein, ignore.case = TRUE) | 
                                 grepl(ground_truth_organisms[2], result.read$Protein, ignore.case = TRUE), 1, 0)
        } else {
            ground_truth <- ifelse(grepl(ground_truth_organisms[1], rownames(result.read), ignore.case = TRUE) | 
                                 grepl(ground_truth_organisms[2], rownames(result.read), ignore.case = TRUE), 1, 0)
        }
    } else {
        # Single organism expected to be differentially expressed
        ground_truth_name <- ground_truth_organisms[1]
        
        if(dea.tool == "ROTS"){
            ground_truth <- ifelse(grepl(ground_truth_name, rownames(result.read$data), ignore.case = TRUE), 1, 0)
            # Fallback for "UPS" vs. "HUMAN" identifier mismatch
            if(ground_truth_name == "UPS"){
                ground_truth_human <- ifelse(grepl("HUMAN", rownames(result.read$data), ignore.case = TRUE), 1, 0)
                if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
                    ground_truth <- ground_truth_human
                }
            }
        } else if(dea.tool == "MSstats"){
            ground_truth <- ifelse(grepl(ground_truth_name, result.read$Protein, ignore.case = TRUE), 1, 0)
            # Fallback for "UPS" vs. "HUMAN" identifier mismatch
            if(ground_truth_name == "UPS"){
                ground_truth_human <- ifelse(grepl("HUMAN", result.read$Protein, ignore.case = TRUE), 1, 0)
                #ground_truth_notEcoli <- ifelse(!grepl("ECOLI", result.read$Protein, ignore.case = TRUE), 1, 0)
                
                if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
                    ground_truth <- ground_truth_human
                }
                

            }
        } else {
            ground_truth <- ifelse(grepl(ground_truth_name, rownames(result.read), ignore.case = TRUE), 1, 0)
            # Fallback for "UPS" vs. "HUMAN" identifier mismatch
            if(ground_truth_name == "UPS"){
                ground_truth_human <- ifelse(grepl("HUMAN", rownames(result.read), ignore.case = TRUE), 1, 0)
                if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
                    ground_truth <- ground_truth_human
                }
            }
        }
    }
    
    # If no positive/negative examples, attempt to remap protein identifiers
    if(length(unique(ground_truth)) != 2){

        
        # Special case: MSstats and DEqMS need Uniprot ID mapping
        if(dea.tool == "MSstats"){
            # Handle file naming: Spectronaut uses "spt" abbreviation in filenames
            analysis_suffix <- if(ANALYSIS_MODE == "Spectronaut") "spt" else ANALYSIS_MODE
            proteins <- read.csv(file.path(DATA_DIR, paste0(exp.name, "_DIA_", analysis_suffix, "_all_proteins.tsv")), sep = "\t")
            ids <- str_replace_all(proteins$Protein, fixed("tr|"), fixed("sp|"))
            extracted_ids <- sapply(str_extract_all(ids, "(?<=sp\\|)[^|]+"), paste, collapse = ";")
            proteins$extracted_ids <- extracted_ids
            result.read <- merge(result.read, proteins, by.x = "Protein", by.y = "extracted_ids")
            result.read$Protein <- result.read$Protein.y
            
            # Recompute ground truth with mapped proteins
            if(length(ground_truth_organisms) == 2){
                ground_truth <- ifelse(grepl(ground_truth_organisms[1], result.read$Protein, ignore.case = TRUE) | 
                                     grepl(ground_truth_organisms[2], result.read$Protein, ignore.case = TRUE), 1, 0)
            } else {
                ground_truth <- ifelse(grepl(ground_truth_organisms[1], result.read$Protein, ignore.case = TRUE), 1, 0)
                # Fallback for UPS vs HUMAN
                if(ground_truth_organisms[1] == "UPS"){
                    ground_truth_human <- ifelse(grepl("HUMAN", result.read$Protein, ignore.case = TRUE), 1, 0)
                    if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
                        ground_truth <- ground_truth_human
                    }
                }
            }

        } else {
            warning(paste("No valid ground truth found for", i.file))
            stop("Critical error in ground truth mapping - stopping execution")
        }
    }

    if(is.matrix(result.read)){
        result.read <- data.frame(result.read, check.rows = F, check.names = F)
    }
    
    # Skip file if result is empty
    if(is.data.frame(result.read)){
        if(nrow(result.read) == 0){
            pass = FALSE
        } else {
            pass = TRUE
        }
    } else {
        pass = TRUE
    }

    # Continue only if valid result
    if(pass){
        # Format result table into a common structure
        if(dea.tool == "ttest"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$pval,
                                    fdr = result.read$adj.pvalue,
                                    logfc = NA,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "Limma"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$P.Value,
                                    fdr = result.read$adj.P.Val,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "SAM"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar.,
                                    fdr = result.read$adj.pvalue,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "ANOVA"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$P.Value,
                                    fdr = result.read$adj.P.Val,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "ROTS"){
            common.df <- data.frame(proteins = rownames(result.read$data),
                                    pvalue = result.read$pvalue,
                                    fdr = result.read$FDR,
                                    logfc = result.read$logfc,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "LimROTS"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read@elementMetadata$pvalue,
                                    fdr = result.read@elementMetadata$FDR,
                                    logfc = result.read@elementMetadata$corrected.logfc,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "DEP"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$p.val,
                                    fdr = result.read$p.adj,
                                    logfc = result.read$diff,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "DEqMS"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$sca.P.Value,
                                    fdr = result.read$sca.adj.pval,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        } else if(dea.tool == "MSstats"){
            common.df <- data.frame(proteins = result.read$Protein,
                                    pvalue = result.read$pvalue,
                                    fdr = result.read$adj.pvalue,
                                    logfc = result.read$log2FC,
                                    ground_truth = ground_truth)
            common.df <- common.df[!is.na(common.df$pvalue),]
        }

        # Define score as 1 - FDR for ranking
        common.df$method_score <- 1 - common.df$fdr

        # Compute confusion matrix components
        TP <- length(which(common.df$ground_truth == 1 & common.df$fdr < fdr_threshold))
        TN <- length(which(common.df$ground_truth == 0 & common.df$fdr > fdr_threshold))
        FP <- length(which(common.df$ground_truth == 0 & common.df$fdr < fdr_threshold))
        FN <- length(which(common.df$ground_truth == 1 & common.df$fdr > fdr_threshold))

        # Convert to numeric to prevent overflow in integer operations
        TP <- as.numeric(TP)
        TN <- as.numeric(TN)
        FP <- as.numeric(FP)
        FN <- as.numeric(FN)

        # Compute partial AUC if ground truth is valid
        if(length(unique(common.df$ground_truth)) >= 2){
            roc_curve <- roc(common.df$ground_truth, common.df$method_score, smooth = F)  
            pauc <- auc(roc_curve, partial.auc = c(1 - fdr_threshold, 1),
                        partial.auc.correct = TRUE)
        } else {
            pauc <- NA
        }

        # Helper function for safe division
        safe_div <- function(numerator, denominator) {
            if (denominator > 0) return(numerator / denominator)
            return(NA)
        }

        metrics <- data.frame(
            Metric = c(
                "TP", "TN", "FP", "FN",
                "Accuracy", 
                "Precision", 
                "Recall (Sensitivity)", 
                "Specificity", 
                "F1_Score", 
                "False_Discovery_Rate", 
                "False_Negative_Rate", 
                "Matthews_Correlation_Coefficient",
                "pauc",
                "G_Mean",
                "Fowlkes_Mallows_Index",
                "Balanced_Accuracy"
            ),
            Value = c(
                TP, TN, FP, FN,
                # Accuracy
                safe_div((TP + TN), (TP + TN + FP + FN)),
                # Precision
                safe_div(TP, (TP + FP)),
                # Recall (Sensitivity)
                safe_div(TP, (TP + FN)),
                # Specificity
                safe_div(TN, (TN + FP)),
                # F1-Score
                safe_div(2 * TP, (2 * TP + FP + FN)),
                # False Discovery Rate
                safe_div(FP, (TP + FP)),
                # False Negative Rate
                safe_div(FN, (TP + FN)),
                # Matthews Correlation Coefficient (MCC)
                ifelse(
                    (TP + FP) * (TP + FN) * (TN + FP) * (TN + FN) > 0,
                    ((TP * TN) - (FP * FN)) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)),
                    NA
                ),
                # pAUC
                pauc,
                # G-Mean
                sqrt(safe_div(TP, (TP + FN)) * safe_div(TN, (TN + FP))),
                # Fowlkes-Mallows Index
                sqrt(safe_div(TP, (TP + FP)) * safe_div(TP, (TP + FN))),
                # Balanced Accuracy
                safe_div(safe_div(TP, (TP + FN)) + safe_div(TN, (TN + FP)), 2)
            )
        )
        
        metrics[17,1] <- "nMCC"
        metrics[17,2] <- (1/2)*(metrics[12,2] + 1)
        
        colnames(metrics)[2] <- i.file
        
        row.names(metrics) <- metrics$Metric
        metrics$Metric <- NULL
        
        metrics.all <- cbind(metrics, metrics.all)
    }
}

# Save results with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_file <- file.path(OUTPUT_DIR, sprintf("metrics_%s_FDR%.2f_%s.csv", 
                                           "Case4", fdr_threshold, timestamp))
write.csv(metrics.all, output_file)

message(sprintf("Analysis completed for %s mode\nResults saved to: %s", 
               ANALYSIS_MODE, output_file))
