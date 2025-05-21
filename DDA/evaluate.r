# Load required libraries
library(stringr)  
library(pROC)    
library(caret) 

# List files in the data folder
data.file <- list.files("data/")

# Load experimental design info
data.info <- read.csv("data/data.info.csv")

# Extract experiment ID from dataset names
data.info$Dataset <- str_split_fixed(data.info$Dataset , fixed("_"), 2)[,1]

# List differential expression result files
results.files <- list.files("Maxquant_results//")

# Set FDR threshold for classification
fdr_threshold <- 0.05

# Initialize dataframe to store performance metrics
metrics.all <- data.frame(row.names = 1:17)

# Loop over each result file
for(i.file in results.files){
    # Extract tool, experiment, and contrast from filename
    get.info <- str_split_fixed(i.file , fixed("_"), 3)
    dea.tool <- get.info[,1]
    exp.name <- get.info[,2]
    Contrast <- get.info[,3]
    
    # Extract condition from contrast
    wide_narrow <- str_split_fixed(Contrast,  fixed("_") , 2)[,1]
    
    # Refine Contrast by removing prefix and suffix
    if(grepl("_" , Contrast, fixed = TRUE)){
        Contrast <- str_split_fixed(Contrast, fixed("_"), 2)[,2]
    }
    Contrast <- str_remove(Contrast, fixed(".rds"))

    # Get the design file name for this experiment
    design <- data.info[data.info$Dataset == exp.name,]
    design <- unique(design$ID)
    design <- which(sapply(data.file, function(x) grepl(design, x, fixed = TRUE) ) == TRUE)
    design <- names(design)
    design <- read.csv(paste0("data/", design), sep = "\t")
    
    # Load DEA results
    result.read <- readRDS(paste0("Maxquant_results//" , i.file))
    ground_truth_name <- colnames(design)[2]
    
    # Determine ground truth for each DEA tool
    if(dea.tool == "ROTS"){
        ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read$data) ) )
        ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read$data)  ) ,1 , 0)
    }else{
        ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read) ) )
        ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read)  ) ,1 , 0)
    }
    
    if(dea.tool == "MSstats"){
        ground_truth_id <- which(grepl(ground_truth_name , result.read$Protein ) )
        ground_truth <- ifelse( grepl(ground_truth_name , result.read$Protein   ) ,1 , 0)
    }
    
    # Fallback for "UPS" vs. "HUMAN" identifier mismatch
    if(length(ground_truth_id) == 0){
      if(tolower(ground_truth_name) == "ups"){
        ground_truth_name <- "HUMAN"
        if(dea.tool == "ROTS"){
          ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read$data) ) )
          ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read$data)  ) ,1 , 0)
        }else{
          ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read) ) )
          ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read)  ) ,1 , 0)
        }
        
        if(dea.tool == "MSstats"){
          ground_truth_id <- which(grepl(ground_truth_name , result.read$Protein ) )
          ground_truth <- ifelse( grepl(ground_truth_name , result.read$Protein   ) ,1 , 0)
        }
      }
    }
    
    # If no positive/negative examples, attempt to remap protein identifiers
    if(length(unique(ground_truth)) != 2){
        
        if(exp.name == "HEof"){
            exp.name_ <- paste0(exp.name , "_" , wide_narrow)
        }else{
            exp.name_ <- exp.name
        }

        # Special case: MSstats and DEqMS need Uniprot ID mapping
        if(dea.tool == "MSstats"){
            proteins <- read.csv( paste0( "../Maxquant/" , exp.name_ , "_DIA_spt_all_proteins.tsv"), sep = "\t" )
            ids <- str_replace_all(proteins$Protein, fixed("tr|"), fixed("sp|"))
            extracted_ids <- sapply(str_extract_all(ids, "(?<=sp\\|)[^|]+"), paste, collapse = ";")
            proteins$extracted_ids <- extracted_ids
            result.read <- merge(result.read , proteins, by.x = "Protein" , by.y = "extracted_ids")
            result.read$Protein <- result.read$Protein.y
            ground_truth_id <- which(grepl(ground_truth_name , result.read$Protein ) )
            ground_truth <- ifelse( grepl(ground_truth_name , result.read$Protein   ) ,1 , 0)

        }else if(dea.tool == "DEqMS"){
            proteins <- read.csv( paste0( "../Maxquant/" , exp.name_ , "_DIA_spt_all_proteins.tsv"), sep = "\t" )
            ids <- str_replace_all(proteins$Protein, fixed("tr|"), fixed("sp|"))
            extracted_ids <- sapply(str_extract_all(ids, "(?<=sp\\|)[^|]+"), paste, collapse = ";")
            proteins$extracted_ids <- extracted_ids
            result.read$Protein <- row.names(result.read)
            result.read <- merge(result.read , proteins, by.x = "Protein" , by.y = "extracted_ids")
            result.read$Protein <- result.read$Protein.y
            ground_truth_id <- which(grepl(ground_truth_name , result.read$Protein ) )
            ground_truth <- ifelse( grepl(ground_truth_name , result.read$Protein   ) ,1 , 0)

        }else{
            stop("No ground_truth")
        }
    }

    # Skip file if result is empty
    if(is.data.frame(result.read)){
        if(nrow(result.read) == 0){
            pass = FALSE
        }else{
            pass = TRUE
        }
    }else{
        pass = TRUE
    }

    # Continue only if valid result
    if(pass){
        # Format result table into a common structure
        if(dea.tool == "ttest"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$pvalue,
                                    fdr = result.read$adj.pvalue,
                                    logfc = result.read$fc.calc,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "Limma"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$P.Value,
                                    fdr = result.read$adj.P.Val,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "SAM"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar.,
                                    fdr = result.read$adj.pvalue,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "ANOVA"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$P.Value,
                                    fdr = result.read$adj.P.Val,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "ROTS"){
            common.df <- data.frame(proteins = rownames(result.read$data),
                                    pvalue = result.read$pvalue,
                                    fdr = result.read$FDR,
                                    logfc = result.read$logfc,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "LimROTS"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read@elementMetadata$pvalue,
                                    fdr = result.read@elementMetadata$FDR,
                                    logfc = result.read@elementMetadata$corrected.logfc,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "DEP"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$p.val,
                                    fdr = result.read$p.adj,
                                    logfc = result.read$diff,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "DEqMS"){
            common.df <- data.frame(proteins = rownames(result.read),
                                    pvalue = result.read$sca.P.Value,
                                    fdr = result.read$sca.adj.pval,
                                    logfc = result.read$logFC,
                                    ground_truth = ground_truth)
        }else if(dea.tool == "MSstats"){
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
        TN <- length(which(common.df$ground_truth == 0 & common.df$fdr > fdr_threshold  ))
        FP <- length(which(common.df$ground_truth == 0 & common.df$fdr < fdr_threshold))
        FN <- length(which(common.df$ground_truth == 1 & common.df$fdr > fdr_threshold  ))

        # Convert to numeric to prevent overflow in integer operations
        TP <- as.numeric(TP)
        TN <- as.numeric(TN)
        FP <- as.numeric(FP)
        FN <- as.numeric(FN)

        # Compute partial AUC if ground truth is valid
        if(length(unique(common.df$ground_truth)) >= 2){
            roc_curve <- roc(common.df$ground_truth, common.df$method_score , smooth = F)  
            pauc <- auc(roc_curve, partial.auc = c(1 - fdr_threshold, 1 ),
                        partial.auc.correct=TRUE)
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
            "Balanced_Accuracy"  # Added Balanced Accuracy
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
            # pAUC (placeholder)
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
    
    metrics.all <- cbind(metrics , metrics.all)
    


    }
}
write.csv(metrics.all , paste0("metrics_Maxquant//" , "FDR_" , fdr_threshold ))

