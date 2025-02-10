library(stringr)
library(pROC)
library(caret)

data.file <- list.files("data/")
data.info <- read.csv("data/data.info.csv")
data.info$Dataset <- str_split_fixed(data.info$Dataset , fixed("_"), 2)[,1]



results.files <- list.files("FragPipe_results/")

true_fold <- function(Contrast, design){
    ground_truth_name <- colnames(design)[2]
    #design <- design[order(design[,ground_truth_name],decreasing = T),]
    Contrast.df <- design[design$condition %in% c(str_split(Contrast, "")[[1]]),]
    Contrast.df <- Contrast.df[match(Contrast.df$condition, c(str_split(Contrast, "")[[1]])),]
    true.fc <- Contrast.df[,ground_truth_name][1] / Contrast.df[,ground_truth_name][2]
    true.fc <- log2(true.fc)
    return(true.fc)
}



fdr_threshold <- 0.05
logfc_threshold <- log2(1)

metrics.all <- data.frame(row.names = 1:17)

for(i.file in results.files){
    get.info <- str_split_fixed(i.file , fixed("_"), 3)
    dea.tool <- get.info[,1]
    exp.name <- get.info[,2]
    Contrast <- get.info[,3]
    
    if(grepl("_" , Contrast, fixed = TRUE)){
        Contrast <- str_split_fixed(Contrast, fixed("_"), 2)[,2]
    }
    Contrast <- str_remove(Contrast, fixed(".rds"))
    design <- data.info[data.info$Dataset == exp.name,]
    design <- unique(design$ID)
    design <- which(sapply(data.file, function(x) grepl(design, x, fixed = TRUE) ) == TRUE)
    design <- names(design)
    design <- read.csv(paste0("data/", design), sep = "\t")
    true.fc <- true_fold(Contrast, design)
    
    result.read <- readRDS(paste0("FragPipe_results/" , i.file))
    ground_truth_name <- colnames(design)[2]
    
    if(dea.tool == "ROTS"){
        ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read$data) ) )
        ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read$data)  ) ,1 , 0)
        if(length(ground_truth_id) == 0){
            ground_truth_id <- which(grepl("HUMAN" , rownames(result.read$data) ) )
            ground_truth <- ifelse( grepl("HUMAN" , rownames(result.read$data) ) ,1 , 0)
            
        }
    }else{
    ground_truth_id <- which(grepl(ground_truth_name , rownames(result.read) ) )
    ground_truth <- ifelse( grepl(ground_truth_name , rownames(result.read)  ) ,1 , 0)
    
    if(length(ground_truth_id) == 0){
        ground_truth_id <- which(grepl("HUMAN" , rownames(result.read) ) )
        ground_truth <- ifelse( grepl("HUMAN" , rownames(result.read)  ) ,1 , 0)
        
    }
    
    }
    
    if(dea.tool == "ttest"){
        
        common.df <- data.frame(proteins = rownames(result.read),
                                pvalue = result.read$pvalue,
                                fdr = result.read$adj.pvalue,
                                logfc = result.read$fc.calc,
                                ground_truth = ground_truth)
        
        if(true.fc < 0){
            common.df$logfc <- -1 * common.df$logfc
        }
        
        
    }else if(dea.tool == "Limma"){
        
        common.df <- data.frame(proteins = rownames(result.read),
                                pvalue = result.read$P.Value,
                                fdr = result.read$adj.P.Val,
                                logfc = result.read$logFC,
                                ground_truth = ground_truth)
        if(true.fc < 0){
            common.df$logfc <- -1 * common.df$logfc
        }        
        
    }else if(dea.tool == "SAM"){
        
        common.df <- data.frame(proteins = rownames(result.read),
                                pvalue = result.read$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar.,
                                fdr = result.read$adj.pvalue,
                                logfc = result.read$logFC,
                                ground_truth = ground_truth)
        
        if(true.fc >= 0){
            common.df$logfc <- -1 * common.df$logfc
        }        
        
    }else if(dea.tool == "ANOVA"){
        
        common.df <- data.frame(proteins = rownames(result.read),
                                pvalue = result.read$P.Value,
                                fdr = result.read$adj.P.Val,
                                logfc = result.read$logFC,
                                ground_truth = ground_truth)
        
        if(true.fc >= 0){
            common.df$logfc <- -1 * common.df$logfc
        }        
    }else if(dea.tool == "ROTS"){
        
        common.df <- data.frame(proteins = rownames(result.read$data),
                                pvalue = result.read$pvalue,
                                fdr = result.read$FDR,
                                logfc = result.read$logfc,
                                ground_truth = ground_truth)
        
        if(true.fc < 0){
            common.df$logfc <- -1 * common.df$logfc
        }        
    }else if(dea.tool == "LimROTS"){
        

        common.df <- data.frame(proteins = rownames(result.read),
                                pvalue = result.read@elementMetadata$pvalue,
                                fdr = result.read@elementMetadata$FDR,
                                logfc = result.read@elementMetadata$corrected.logfc,
                                ground_truth = ground_truth)
        
        if(true.fc < 0){
            common.df$logfc <- -1 * common.df$logfc
        }        
        
    }
    
    
    
    #### No need to change #### 
    FC.ground_truth <- str_split_fixed(i.file, fixed("_") , 2)[,2]
    FC.ground_truth <- paste0("ttest_",FC.ground_truth)
    FC.ground_truth <- readRDS(paste0("FragPipe_results/" , FC.ground_truth))
    
    
    if(true.fc < 0){
        FC.ground_truth$fc.calc <- -1 * FC.ground_truth$fc.calc
    }       
    
    FC.ground_truth <- FC.ground_truth[common.df$proteins,]
    common.df$FC.ground_truth <- FC.ground_truth$fc.calc
    #####
    common.df$ground_truth  <- ifelse(common.df$ground_truth == 1 & 
                                          common.df$FC.ground_truth >= 0, 1, 0)
    
    
    common.df$method_score <- 1 - common.df$fdr
    ground_truth_names <- common.df$proteins[ground_truth_id]
    
    TP <- length(which(common.df$ground_truth == 1 & common.df$fdr < fdr_threshold))
    
    TN <- length(which(common.df$ground_truth == 0 & common.df$fdr > fdr_threshold  ))
    
    FP <- length(which(common.df$ground_truth == 0 & common.df$fdr < fdr_threshold))
    
    FN <- length(which(common.df$ground_truth == 1 & common.df$fdr > fdr_threshold  ))
    
    # Convert all inputs to numeric to avoid integer overflow
    TP <- as.numeric(TP)
    TN <- as.numeric(TN)
    FP <- as.numeric(FP)
    FN <- as.numeric(FN)
    
    
    
    if(length(unique(common.df$ground_truth)) >= 2){
        
        roc_curve <- roc(common.df$ground_truth, common.df$method_score , smooth = F)  
        
        pauc <- auc(roc_curve, partial.auc = c(1 - fdr_threshold, 1 ),
                    partial.auc.correct=TRUE)
        
        saveRDS(roc_curve, paste0("roc_curve_objects_FragPipe/" , i.file))
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

write.csv(metrics.all , paste0("metrics_FragPipe/" , "FDR_" , fdr_threshold ))
 