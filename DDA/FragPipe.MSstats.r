library(stringr)
library(SummarizedExperiment)
library(LimROTS)
library(ROTS)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(limma)
library(samr)
library(MSstats)
library(DEqMS)
library(matrixStats)

# List experiment folders (excluding specific ones)
all.exp  <- list.files("../../")
all.exp <- all.exp[!all.exp %in% c("FragPipe" , "case study 1", "Maxquant")]

# Set random seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

#### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### 
#### Loop over experiments
#### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### 

for(i.exp in all.exp){
    print(i.exp)
    
    # Read FragPipe output file
    exp.mar <- read.csv(paste0("../../", i.exp, "/FragPipe/TOP0/combined_protein.tsv"), sep = "\t")
    
    # Read design file
    design <- read.csv(paste0("../../FragPipe/", i.exp , "_FragPipe_design.tsv"), sep = "\t")
    
    # Read MSstats formatted evidence file
    evidence <- read.csv(paste0("../../", i.exp, "/FragPipe/TOP0/MSstats.csv"))
    
    # Read MSstats annotation
    annot <- read.csv(paste0("../../FragPipe/", i.exp , "_FragPipe_design_msstats.tsv"), sep = "\t")
    
    # Generate all valid pairwise condition contrasts
    Contrasts <- combn(design$condition, 2, simplify = F)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
    
    # Format protein data
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    
    # Clean sample names in design and data
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Convert FragPipe evidence and annotation to MSstats format
    raw <- FragPipetoMSstatsFormat(
        annotation = annot, 
        input = evidence)    
    
    # Save raw MSstats-formatted data
    saveRDS(raw, paste0("rawMSstats/", "RawMSstats_" ,i.exp))
    
    # Preprocess data for differential testing
    QuantData <- dataProcess(raw)
    
    # Loop over each contrast
    for(i.Contrasts in Contrasts){
        print(i.Contrasts)
        
        ########################  ########################    ######################## 
        ######################## MSstats
        ########################  ########################    ######################## 
        
        # Extract contrast conditions
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        # Build contrast matrix
        comparison_r <- levels(QuantData$ProteinLevelData$GROUP)
        comparison_r <- ifelse(comparison_r == i.Contrasts1 , -1,
                               ifelse(comparison_r == i.Contrasts2 , 1 , 0))
        
        comparison <- matrix(comparison_r, nrow = 1)
        colnames(comparison) <- levels(QuantData$ProteinLevelData$GROUP)
        row.names(comparison) <- paste0(i.Contrasts1, "-" , i.Contrasts2)
        print(comparison)
        
        # Run group comparison with MSstats
        testResultOneComparison <- groupComparison(contrast.matrix = comparison, data = QuantData)
        
        # Save result table
        saveRDS(testResultOneComparison$ComparisonResult, paste0("FragPipe_results/", "MSstats_" , i.exp , "_" , i.Contrasts, ".rds"))
        
        # Clean up memory
        remove(testResultOneComparison)
        gc()
    }
    
    # Remove objects to free memory
    remove(raw)
    remove(evidence)
    remove(QuantData)
    gc()
}
