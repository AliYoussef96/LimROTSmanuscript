# Configuration
#------------------------------------------------------------------------------
# Set analysis mode: "DIANN" or "Spectronaut"
ANALYSIS_MODE <- "DIANN"  # Change this to switch between DIANN and Spectronaut

# Set paths based on analysis mode
RESULTS_DIR <- if(ANALYSIS_MODE == "DIANN") "DIANN_results" else "Spectronaut_results"
DATA_DIR <- if(ANALYSIS_MODE == "DIANN") "../DIANN" else "../Spectronaut"
OUTPUT_DIR <- if(ANALYSIS_MODE == "DIANN") "metrics_DIANN" else "metrics_Spectronaut"

# Create output directory if it doesn't exist
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#------------------------------------------------------------------------------
# Load required libraries
suppressPackageStartupMessages({
    library(stringr)
    library(pROC)
    library(caret)
    library(PRROC)
    library(purrr)
    library(dplyr)
    library(ggplot2)
    library(ggsci)
    library(gridExtra)  # For combining plots
    library(grid)       # For plot annotations
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

# Set thresholds
fdr_threshold <- 0.05

# Initialize dataframe to store all results for ROC/PR curves
metrics.all <- data.frame()

# Loop over each result file
for(i.file in results.files){
    # Extract tool, experiment, and contrast from filename
    
    get.info <- str_split_fixed(i.file, fixed("_"), 3)
    dea.tool <- get.info[,1]
    exp.name <- get.info[,2]
    Contrast <- get.info[,3]
    Contrast <- str_remove(Contrast , ".rds")
    

    # Get ground truth organisms for this experiment
    ground_truth_organisms <- get_ground_truth_DIA(exp.name)
    
    if(is.null(ground_truth_organisms)) {
        warning(paste("Skipping", i.file, "- no ground truth found for", exp.name))
        stop("no ground truth found for")
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
            #ground_truth_notEcoli <- ifelse(!grepl("ECOLI", result.read$Protein, ignore.case = TRUE), 1, 0)
            
            # Fallback for "UPS" vs. "HUMAN" identifier mismatch
            if(ground_truth_name == "UPS"){
                ground_truth_human <- ifelse(grepl("HUMAN", result.read$Protein, ignore.case = TRUE), 1, 0)
                if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
                    ground_truth <- ground_truth_human
                }

            }
        } else {
            ground_truth <- ifelse(grepl(ground_truth_name, rownames(result.read), ignore.case = TRUE), 1, 0)
            # Fallback for "UPS" vs. "HUMAN" identifier mismatch
            # if(ground_truth_name == "UPS"){
            #     ground_truth_human <- ifelse(grepl("HUMAN", rownames(result.read), ignore.case = TRUE), 1, 0)
            #     if(sum(ground_truth_human) >= sum(ground_truth) || sum(ground_truth) == 0) {
            #         ground_truth <- ground_truth_human
            #     }
            # }
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
    
    # Skip file if result is empty or no valid ground truth
    if(is.data.frame(result.read)){
        if(nrow(result.read) == 0 || length(unique(ground_truth)) < 2){
            stop("Critical error in ground truth mapping - stopping execution")
        }
    }
    
    # Format result table into a common structure for ROC/PR analysis
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
    
    # Skip if no valid data
    if(nrow(common.df) == 0) stop(i.file)
    
    # Create score for ranking (1 - p-value for better discrimination)
    common.df$score <- 1 - common.df$pvalue
    
    # Add metadata
    common.df$dea.tool <- dea.tool
    common.df$exp <- exp.name
    common.df$Contrast <- Contrast
    

    metrics.all <- rbind(metrics.all, common.df)
    # Only keep if we have valid ground truth
    # if(length(unique(common.df$ground_truth)) >= 2){
    #     metrics.all <- rbind(metrics.all, common.df)
    # }
}

#------------------------------------------------------------------------------

metrics.all <- metrics.all[metrics.all$dea.tool %in% c("DEP", "Limma", "LimROTS", "ANOVA"),]


# Generate PR Curves
if(nrow(metrics.all) > 0){
    # Calculate PR curves for each method
    
    pr_data <- metrics.all %>%
      group_by(dea.tool) %>%
      group_modify(~{
        fg <- .x$score[.x$ground_truth == 1]
        bg <- .x$score[.x$ground_truth == 0]
        pr <- pr.curve(scores.class0 = fg, scores.class1 = bg, curve = TRUE)
        tibble(
          recall = pr$curve[,1],
          precision = pr$curve[,2],
          dea_tool = .y$dea.tool,
          aucpr = pr$auc.integral
        )
      })

    
    
    # Add method labels with AUC values
    pr_data <- pr_data %>%
        group_by(dea.tool) %>%
        mutate(`Method(AUC)` = paste0(dea.tool, " (", sprintf("%.2f", dplyr::first(aucpr)), ")")) %>%
        ungroup()
    lancet_colors <- pal_lancet()(9)  # Get lancet palette colors
    

        tool_colors <- c(
            lancet_colors[1],      
             lancet_colors[2],    
             lancet_colors[5],  
             lancet_colors[4]     
        )
        
        names(tool_colors) <- unique(pr_data$`Method(AUC)`)

    # Generate PR curve plot
    if(nrow(pr_data) > 0){
        p_pr <- ggplot(pr_data, aes(x = recall, y = precision, group = `Method(AUC)`, color = `Method(AUC)`)) +
          geom_line(size = 1.2) +
          theme_bw() +
          scale_color_lancet() +
          labs(x = "Recall", y = "Precision", title = "A) Precision-Recall Curve") +
          theme(
            text = element_text(color = "black", face = "bold", size = 14),
            axis.title = element_text(color = "black", face = "bold", size = 14),
            axis.text = element_text(color = "black", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(color = "black", face = "bold", size = 11),
            legend.position = "right",
            legend.box = "horizontal",
            plot.title = element_text(color = "black", face = "bold", hjust = 0, size = 16),
            panel.grid.minor = element_blank(),
            plot.margin = margin(10, 10, 10, 10)
          ) +
          guides(color = guide_legend(title = NULL, ncol = 1, byrow = TRUE)) +
            scale_fill_manual(values = tool_colors) +             
            scale_color_manual(values = tool_colors)
        
        # Print AUC summary
        auc_summary_pr <- pr_data %>%
          group_by(dea.tool) %>%
          summarise(mean_aucpr = mean(aucpr, na.rm = TRUE)) %>%
          arrange(desc(mean_aucpr))
        
        print("AUC-PR Summary:")
        print(auc_summary_pr)
    }
    
    #--------------------------------------------------------------------------
    # Generate ROC Curves
    roc_data <- metrics.all %>%
      group_by(dea.tool) %>%
      group_modify(~{
        if(length(unique(.x$ground_truth)) < 2) return(tibble())
        
        roc_obj <- tryCatch({
          roc(response = .x$ground_truth, predictor = .x$score, quiet = TRUE)
        }, error = function(e) return(NULL))
        
        if(is.null(roc_obj) || is.atomic(roc_obj)) return(tibble())
        
        tibble(
          fpr = 1 - rev(roc_obj$specificities),
          tpr = rev(roc_obj$sensitivities),
          auc = as.numeric(auc(roc_obj)),
          dea_tool = .y$dea.tool
        )
      })
    
    # Add method labels with AUC values for ROC
    roc_data <- roc_data %>%
      group_by(dea.tool) %>%
      mutate(`Method(AUC)` = paste0(dea.tool, " (", sprintf("%.2f", dplyr::first(auc)), ")")) %>%
      ungroup()
    
    tool_colors <- c(
        lancet_colors[1],      
        lancet_colors[2],    
        lancet_colors[5],  
        lancet_colors[4]     
    )
    
    names(tool_colors) <- unique(roc_data$`Method(AUC)`)
    
    
    # Generate ROC curve plot
    if(nrow(roc_data) > 0){
        p_roc <- ggplot(roc_data, aes(x = fpr, y = tpr, group = `Method(AUC)`, color = `Method(AUC)`)) +
          geom_line(size = 1.2) +
          geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50", alpha = 0.8) +
          theme_bw() +
          scale_color_lancet() +
          labs(x = "False Positive Rate", y = "True Positive Rate", title = "B) ROC Curve") +
          theme(
            text = element_text(color = "black", face = "bold", size = 14),
            axis.title = element_text(color = "black", face = "bold", size = 14),
            axis.text = element_text(color = "black", face = "bold", size = 12),
            legend.title = element_blank(),
            legend.text = element_text(color = "black", face = "bold", size = 11),
            legend.position = "right",
            legend.box = "horizontal",
            plot.title = element_text(color = "black", face = "bold", hjust = 0, size = 16),
            panel.grid.minor = element_blank(),
            plot.margin = margin(10, 10, 10, 10)
          ) +
          guides(color = guide_legend(title = NULL, ncol = 1, byrow = TRUE)) +
            scale_fill_manual(values = tool_colors) +             
            scale_color_manual(values = tool_colors)              
          
        
        # Print AUC summary
        auc_summary_roc <- roc_data %>%
          group_by(dea.tool) %>%
          summarise(mean_auc = mean(auc, na.rm = TRUE)) %>%
          arrange(desc(mean_auc))
        
        print("AUC-ROC Summary:")
        print(auc_summary_roc)
        
        #----------------------------------------------------------------------
        # Create combined publication-ready figure
        if(exists("p_pr") && exists("p_roc")){
            

            # Save individual plots as well for flexibility
            ggsave(file.path(OUTPUT_DIR, "PR_curve_individual.png"), p_pr, width = 8, height = 5, dpi = 300)
            ggsave(file.path(OUTPUT_DIR, "ROC_curve_individual.png"), p_roc, width = 10, height = 8, dpi = 300)
        }
    }
    
} else {
    warning("No valid data found for ROC/PR curve generation")
}

message(sprintf("Analysis completed for %s mode", ANALYSIS_MODE))
