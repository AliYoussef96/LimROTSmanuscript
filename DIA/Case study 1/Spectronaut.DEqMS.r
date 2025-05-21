# Load required libraries
library(stringr)
library(SummarizedExperiment)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(limma)
library(MSstats)
library(DEqMS)
library(data.table)
library(tidyr)

# List all experiment directories except the ones excluded
all.exp  <- list.files("../")
all.exp <- all.exp[!all.exp %in% c("DIANN" , "case study 1", "Spectronaut")]

# Set seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

# Iterate over each experiment
for(i.exp in all.exp){
  print(i.exp)
  
  # Find the Spectronaut report file
  evidence <- list.files(paste0("../", i.exp, "/Spectronaut/"  ))
  evidence <- evidence[grepl("Report.tsv" , evidence)]
  
  # Read the evidence file (Spectronaut output)
  evidence <- fread(paste0("../", i.exp, "/Spectronaut/" , evidence  ), sep = "\t", nThread = 15)
  
  # Load the design file (condition and replicate info)
  design <- read.csv(paste0("../Spectronaut/", i.exp , "_spt_design.tsv"), sep = "\t")
  
  # Load the log2-transformed data file
  exp.mar <- read.csv(paste0("../Spectronaut/", i.exp , "_spt_dlfq.tsv"), sep = "\t")
  
  # Load MSstats-compatible annotation file
  annot <- read.csv(paste0("../Spectronaut/", i.exp , "_spt_design_msstats.tsv"), sep = "\t")

  # Create all unique pairwise condition contrasts (excluding comparisons within the same condition)
  Contrasts <- combn(design$condition, 2, simplify = F)
  Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
  Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
  
  # Prepare data matrix by setting rownames and cleaning column names
  row.names(exp.mar) <- exp.mar$Protein
  exp.mar$Protein <- NULL
  exp.mar$Organism <- NULL
  design$sample_name <- str_remove(design$sample_name, fixed("_"))
  colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
  
  # Iterate through each contrast
  for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    # Extract condition identifiers from contrast string
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    # Filter design for current contrast conditions
    design.temp <- design[design$condition %in% c(i.Contrasts1,i.Contrasts2 ),]

    # Filter evidence table for relevant conditions
    evidence_temp <- evidence[evidence$R.Condition %in% design.temp$condition,]
    
    # Extract relevant columns: condition, replicate, protein name, and quantity
    df.LFQ <- evidence_temp[,c("R.Condition" , "R.Replicate" , "PG.ProteinNames", "PG.Quantity")]

    # If PG.ProteinNames is missing, use PG.ProteinGroups instead
    if(any(df.LFQ$PG.ProteinNames == "")){
      df.LFQ <- evidence_temp[,c("R.Condition" , "R.Replicate" , "PG.ProteinGroups", "PG.Quantity")]
      colnames(df.LFQ)[3] <- "PG.ProteinNames"
      evidence_temp <- evidence_temp[!duplicated(evidence_temp[,c("PG.ProteinGroups","PEP.StrippedSequence")])]
      pep.count.table <- evidence_temp %>% group_by(PG.ProteinGroups) %>% summarise(count = n())
    } else {
      evidence_temp <- evidence_temp[!duplicated(evidence_temp[,c("PG.ProteinNames","PEP.StrippedSequence")])]
      pep.count.table <- evidence_temp %>% group_by(PG.ProteinNames) %>% summarise(count = n())
    }
    
    # Create unique sample identifiers
    df.LFQ$sample <-  paste0(df.LFQ$R.Condition, df.LFQ$R.Replicate)
    df.LFQ <- df.LFQ[,-c(1,2)]
    
    # Remove duplicate rows by protein and sample
    df.LFQ <- df.LFQ[!duplicated(df.LFQ[,c(1,3)]),]
    
    # Reshape data to wide format: proteins x samples
    df.LFQ <- df.LFQ %>%
      pivot_wider(names_from = PG.ProteinNames, values_from = PG.Quantity)
    df.LFQ <- data.frame(df.LFQ, check.rows = F, check.names = F)
    
    # Set sample names as rownames, transpose to protein x sample matrix
    row.names(df.LFQ) <- df.LFQ$sample
    df.LFQ$sample <- NULL
    df.LFQ <- data.frame(t(df.LFQ), check.rows = F, check.names = F)
    
    # Reconstruct sample names and order
    design.temp$sample_name <- paste0(design.temp$condition, design.temp$replicate)
    design.temp <- design.temp[design.temp$sample_name %in% colnames(df.LFQ),]
    design.temp <- design.temp[order(design.temp$condition),]
    df.LFQ = df.LFQ[,design.temp$sample_name]

    # Prepare peptide count table
    pep.count.table <- data.frame(pep.count.table)
    row.names(pep.count.table) <- pep.count.table[,1]
    
    # Add 1 to peptide counts to avoid 0 in downstream DEqMS step
    pep.count.table$count = pep.count.table$count + 1
    
    # Log2 transform the protein matrix
    protein.matrix = log2(as.matrix(df.LFQ) + 1 )
    
    # Impute missing values using sequential imputation
    protein.matrix = impSeq(protein.matrix)
    
    # Fit linear model without intercept using limma
    design_model = model.matrix(~0+condition, data = design.temp)
    fit1 = lmFit(protein.matrix,design = design_model)
    
    # Create and apply contrast
    cont <- makeContrasts(paste0("condition" , i.Contrasts1, "-condition" , i.Contrasts2), levels = design_model)
    fit2 = contrasts.fit(fit1,contrasts = cont)
    fit3 <- eBayes(fit2)
    
    # Add peptide count info to fit3 object
    fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
    fit3$count[is.na(fit3$count)] <- 1
    
    # Check for invalid peptide counts
    if(min(fit3$count) < 1){
      stop("fit3$count less than1")
    }
    
    # Apply DEqMS correction using spectral count info
    fit4 = spectraCounteBayes(fit3)
    
    # Generate DEqMS results table
    DEqMS.results = outputResult(fit4,coef_col = 1)
    
    # Save result to file
    saveRDS(DEqMS.results, paste0("Spectronaut_results/", "DEqMS_" , i.exp , "_" , i.Contrasts, ".rds"))
    
    # Free memory
    remove(fit3)
    remove(fit4)
    remove(evidence_temp)
    gc()    
  }
  
  # Free memory for current experiment
  remove(evidence)
  gc()
}
