# Load required libraries
library(stringr)
library(MSstats)
library(DEqMS)
library(data.table)

# List all experiments in the parent directory, excluding specific folders
all.exp  <- list.files("../")
all.exp <- all.exp[!all.exp %in% c("DIANN" , "Case 1 DIANN", "Spectronaut" , "case 1 spectronaut")]

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

# Loop over each experiment
for(i.exp in all.exp){
  print(i.exp)
  
  # Locate DIANN report file within the experiment directory
  exp.mar <- list.files(paste0("../", i.exp, "/DIANN/"  ))
  exp.mar <- exp.mar[grepl("report.tsv" , exp.mar)]
  
  # Read DIANN output table using data.table for faster reading
  exp.mar <- fread(paste0("../", i.exp, "/DIANN/" , exp.mar  ), sep = "\t", nThread = 15)
  
  # Read design file (sample conditions, etc.)
  design <- read.csv(paste0("../DIANN/", i.exp , "_DIANN_design.tsv"), sep = "\t")
  
  # Read MSstats-specific annotation file
  annot <- read.csv(paste0("../DIANN/", i.exp , "_DIANN_design_msstats.tsv"), sep = "\t")
  
  # Create unique pairwise contrasts between conditions
  Contrasts <- combn(design$condition, 2, simplify = F)
  Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
  Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
  
  # Clean and align data
  row.names(exp.mar) <- exp.mar$Protein
  exp.mar$Protein <- NULL
  exp.mar$Organism <- NULL
  design$sample_name <- str_remove(design$sample_name, fixed("_"))
  colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
  
  # Convert DIANN format to MSstats format
  raw <- DIANNtoMSstatsFormat(input=exp.mar, 
                             annotation=annot)    
  
  # Save raw MSstats formatted data
  saveRDS(raw, paste0("rawMSstats/", "RawMSstats_" ,i.exp))
  
  # Process data using MSstats (normalization, summarization, etc.)
  QuantData <- dataProcess(raw, numberOfCores = 20)
  
  # Save processed data
  saveRDS(QuantData, paste0("rawMSstats/", "QuantData_" ,i.exp))
  
  # Loop through each contrast
  for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    ########################  ########################    ######################## 
    ######################## MSstats
    ########################  ########################    ######################## 
    
    # Extract group labels from contrast string
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    # Construct contrast matrix for MSstats comparison
    comparison_r <- levels(QuantData$ProteinLevelData$GROUP )
    comparison_r <- ifelse(comparison_r == i.Contrasts1 , -1,
                           ifelse(comparison_r == i.Contrasts2 , 1 , 0))
    
    comparison <- matrix(comparison_r, nrow=1)
    colnames(comparison) <- levels(QuantData$ProteinLevelData$GROUP )
    row.names(comparison) <- paste0(i.Contrasts1, "-" , i.Contrasts2)
    print(comparison)
    
    # Run group comparison with MSstats
    testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData,
                                               numberOfCores = 20, verbose = FALSE,
                                               save_fitted_models = FALSE)
    
    # Save comparison result
    saveRDS(testResultOneComparison$ComparisonResult, paste0("DIANN_results/", "MSstats_" , i.exp , "_" , i.Contrasts, ".rds"))
    
    # Clean up memory
    remove(testResultOneComparison)
    gc()
    
  }
  
  # Clean up memory for next iteration
  remove(raw)
  remove(exp.mar)
  remove(QuantData)
  gc()
}
