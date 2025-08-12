# Load required libraries
library(stringr)    
library(MSstats)   
library(DEqMS)     

# Set experiment identifier
i.exp <- "HEof"

# Read expression and annotation data for the "Narrow" condition
exp.mar1 <- read.csv("../HEof_n600_DIA/DIANN/report.tsv", sep = "\t")
annot.mar1 <- read.csv("../DIANN/HEof_n600_DIA_DIANN_design_msstats.tsv", sep = "\t")
annot.mar1$Condition <- paste0(annot.mar1$Condition , "_Narrow")      # Append "_Narrow" to condition labels
exp.mar1$R.Condition <- paste0(exp.mar1$R.Condition , "_Narrow")      # Append "_Narrow" to raw condition labels

# Read expression and annotation data for the "Wide" condition
exp.mar2 <- read.csv("../HEof_w600_DIA/DIANN/report.tsv", sep = "\t")
annot.mar2 <- read.csv("../DIANN/HEof_w600_DIA_DIANN_design_msstats.tsv", sep = "\t")
annot.mar2$Condition <- paste0(annot.mar2$Condition , "_Wide")        # Append "_Wide" to condition labels
exp.mar2$R.Condition <- paste0(exp.mar2$R.Condition , "_Wide")        # Append "_Wide" to raw condition labels

# Combine "Narrow" and "Wide" datasets
exp <- rbind(exp.mar1, exp.mar2)
annot <- rbind(annot.mar1, annot.mar2)

# Use Protein.Ids column as Protein.Names for MSstats compatibility
exp$Protein.Names <- exp$Protein.Ids

# Convert DIA-NN format to MSstats format
raw <- DIANNtoMSstatsFormat(
    annotation =annot,
    input=exp)    

# Save preprocessed data
saveRDS(raw, paste0("rawMSstats/", "RawMSstats_" ,i.exp))

# Run data processing (normalization, summarization)
QuantData <- dataProcess(raw,  numberOfCores = 20)

# Save the processed data
saveRDS(QuantData, paste0("rawMSstats/", "QuantDataMSstats_" , i.exp))

# Clean up memory
remove(raw)
remove(exp.mar1)
remove(exp.mar2)
gc()  # Trigger garbage collection

###########
# Generate all unique pairwise contrasts between condition labels A to H
# Exclude duplicates and self-contrasts
Contrasts <- combn(c("A", "B" , "C" , "D" , "E" , "F", "G" , "H"), 2, simplify = F)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]

# Set random seed for reproducibility
set.seed(123456, sample.kind = "default" , kind = "default")

# Loop through each pairwise contrast
for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    # Extract the two condition labels
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    # Create contrast labels for Narrow and Wide conditions
    i.Contrasts1.n <- paste0(i.Contrasts1 , "_Narrow")
    i.Contrasts1.w <- paste0(i.Contrasts1 , "_Wide")
    
    i.Contrasts2.n <- paste0(i.Contrasts2 , "_Narrow")
    i.Contrasts2.w <- paste0(i.Contrasts2 , "_Wide")
    
    # Generate contrast vector for MSstats
    comparison_r <- levels(QuantData$ProteinLevelData$GROUP )
    comparison_r <- ifelse(comparison_r == i.Contrasts1.n , 0.5,
                           ifelse(comparison_r == i.Contrasts2.n , -0.5 , 
                                  ifelse(comparison_r == i.Contrasts1.w , 0.5 , 
                                         ifelse(comparison_r == i.Contrasts2.w, -0.5, 0))))
    
    comparison <- matrix(comparison_r, nrow=1)
    colnames(comparison) <- levels(QuantData$ProteinLevelData$GROUP )
    row.names(comparison) <- paste0(i.Contrasts1, "-" , i.Contrasts2)
    print(comparison)
    
    # Perform group comparison using MSstats
    testResultOneComparison <- groupComparison(contrast.matrix=comparison, data=QuantData , numberOfCores = 20)
    
    # Save results
    saveRDS(testResultOneComparison$ComparisonResult, paste0("results_case2_DIAN/", "MSstats_" , i.exp , "_" , i.Contrasts, ".rds"))
}
