# Load required libraries
library(stringr)
library(MSstats)
library(DEqMS)
library(data.table)

# List all experiment folders except specific ones
all.exp  <- list.files("../")
all.exp <- all.exp[!all.exp %in% c("DIANN" , "case study 1", "Spectronaut")]

# Set seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

# Loop over each experiment
for(i.exp in all.exp){
  print(i.exp)
  
  # List Spectronaut report files in the current experiment folder
  exp.mar <- list.files(paste0("../", i.exp, "/Spectronaut/"  ))
  exp.mar <- exp.mar[grepl("Report.tsv" , exp.mar)]  # Keep only Report.tsv file
  
  # Read Spectronaut report data using fread for speed
  exp.mar <- fread(paste0("../", i.exp, "/Spectronaut/" , exp.mar  ), sep = "\t", nThread = 15)
  
  # Read the experimental design
  design <- read.csv(paste0("../Spectronaut/", i.exp , "_spt_design.tsv"), sep = "\t")
  
  # (Optional) Read MaxQuant evidence file - currently commented out
  # evidence <- read.csv(paste0("../", i.exp, "/Spectronaut/evidence.txt"  ), sep = "\t")
  
  # Read MSstats annotation file
  annot <- read.csv(paste0("../Spectronaut/", i.exp , "_spt_design_msstats.tsv"), sep = "\t")
  
  # Create all pairwise contrasts of conditions
  Contrasts <- combn(design$condition, 2, simplify = F)
  Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))  # Collapse to string to avoid duplicates
  Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]  # Remove comparisons with same group (e.g., "AA")

  # Set Protein as row names, and remove non-numeric columns
  row.names(exp.mar) <- exp.mar$Protein
  exp.mar$Protein <- NULL
  exp.mar$Organism <- NULL
  
  # Clean underscores from sample names in design and matrix column names
  design$sample_name <- str_remove(design$sample_name, fixed("_"))
  colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
  
  # Optionally filter columns to match design sample names (currently commented out)
  # exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]

  # Convert Spectronaut format to MSstats format
  raw <- SpectronauttoMSstatsFormat(input=exp.mar, annotation=annot)    
  
  # Save the raw formatted data
  saveRDS(raw, paste0("rawMSstats/", "RawMSstats_" ,i.exp))
  
  # Process the data with MSstats (normalization, summarization, etc.)
  QuantData <- dataProcess(raw, numberOfCores = 20)
  
  # Save the processed data
  saveRDS(QuantData, paste0("rawMSstats/", "QuantData_" ,i.exp))

  # Loop through each contrast to perform group comparison
  for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    # Extract group names from contrast string (assumes 1 char per group)
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    # Create contrast vector: -1 for group1, +1 for group2, 0 otherwise
    comparison_r <- levels(QuantData$ProteinLevelData$GROUP )
    comparison_r <- ifelse(comparison_r == i.Contrasts1 , -1,
                           ifelse(comparison_r == i.Contrasts2 , 1 , 0))
    
    # Create contrast matrix
    comparison <- matrix(comparison_r, nrow=1)
    colnames(comparison) <- levels(QuantData$ProteinLevelData$GROUP )
    row.names(comparison) <- paste0(i.Contrasts1, "-" , i.Contrasts2)
    print(comparison)
    
    # Run MSstats group comparison
    testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData,
                                               numberOfCores = 20, verbose = FALSE,
                                               save_fitted_models = FALSE)
    
    # (Optional) Filter NA p-values - currently commented out
    # testResultOneComparison_df <- testResultOneComparison$ComparisonResult
    # testResultOneComparison_df <- testResultOneComparison_df[!is.na(testResultOneComparison_df$pvalue),]
    
    # Save result of this contrast
    saveRDS(testResultOneComparison$ComparisonResult, paste0("Spectronaut_results/", "MSstats_" , i.exp , "_" , i.Contrasts, ".rds"))
    
    # Clean memory
    remove(testResultOneComparison)
    gc()
  }
  
  # Clean memory
  remove(raw)
  remove(exp.mar)
  remove(QuantData)
  gc()
}
