# Load required libraries for data manipulation, statistical analysis, and parallel processing
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

# List all experiment directories, excluding specified ones
all.exp  <- list.files("../")
all.exp <- all.exp[!all.exp %in% c("Case 1 DIANN" , "case 1 spectronaut", "Spectronaut" , "DIANN")]

# Set seed for reproducibility
set.seed(1597, sample.kind = "default" , kind = "default")

# Iterate over each experiment
for(i.exp in all.exp){
  print(i.exp)
  
  # Identify and read the DIANN report file
  evidence <- list.files(paste0("../", i.exp, "/DIANN/"  ))
  evidence <- evidence[grepl("report.tsv" , evidence)]
  evidence <- fread(paste0("../", i.exp, "/DIANN/" , evidence  ), sep = "\t", nThread = 15)
  
  # Read design and annotation files
  design <- read.csv(paste0("../DIANN/", i.exp , "_DIANN_design.tsv"), sep = "\t")
  exp.mar <- read.csv(paste0("../DIANN/", i.exp , "_DIANN_dlfq.tsv"), sep = "\t")
  annot <- read.csv(paste0("../DIANN/", i.exp , "_DIANN_design_msstats.tsv"), sep = "\t")
  
  # Generate unique pairwise contrasts between conditions
  Contrasts <- combn(design$condition, 2, simplify = F)
  Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
  Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
  
  # Prepare expression matrix
  row.names(exp.mar) <- exp.mar$Protein
  exp.mar$Protein <- NULL
  exp.mar$Organism <- NULL
  design$sample_name <- str_remove(design$sample_name, fixed("_"))
  colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
  
  # Iterate over each contrast
  for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    # Extract group labels from contrast string
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    # Subset design for current contrast
    design.temp <- design[design$condition %in% c(i.Contrasts1,i.Contrasts2 ),]
    
    # Process file names to match 'Run' identifiers in evidence
    design.temp$File.Name <- basename( gsub("\\\\", "/", design.temp$file) ) 
    design.temp$File.Name  <- str_split_fixed(design.temp$File.Name , fixed("."),2)[,1]
    
    # Subset evidence data for current contrast
    evidence_temp <- evidence[evidence$Run %in% design.temp$File.Name,]
    
    # Extract relevant columns for LFQ analysis
    df.LFQ <- evidence_temp[,c("Run" , "Protein.Group" ,
                               "PG.MaxLFQ", "Stripped.Sequence")]
    colnames(df.LFQ)[2] <- "PG.ProteinNames"
    
    # Remove duplicate entries
    evidence_temp <- evidence_temp[!duplicated(evidence_temp[,c("Protein.Group","Stripped.Sequence")])]
    
    # Generate peptide count table
    pep.count.table <- evidence_temp %>% group_by(Protein.Group) %>% summarise(count = n())
    
    # Remove duplicate protein-sample combinations
    df.LFQ <- df.LFQ[!duplicated(df.LFQ[,c("Run","PG.ProteinNames")]),]
    df.LFQ <- df.LFQ[,-4]
    
    # Reshape data to wide format
    df.LFQ <- df.LFQ %>%
      pivot_wider(names_from = PG.ProteinNames, values_from = PG.MaxLFQ)
    df.LFQ <- data.frame(df.LFQ, check.rows = F, check.names = F)
    colnames(df.LFQ)[1] <- "sample"
    row.names(df.LFQ) <- df.LFQ$sample
    df.LFQ$sample <- NULL
    df.LFQ <- data.frame(t(df.LFQ), check.rows = F, check.names = F)
    
    # Proceed if sufficient samples are available
    if(ncol(df.LFQ) >= 5){
      # Align columns with design
      design.temp <- design.temp[order(design.temp$condition),]
      df.LFQ = df.LFQ[,design.temp$File.Name]
      
      # Prepare peptide count table
      pep.count.table <- data.frame(pep.count.table)
      row.names(pep.count.table) <- pep.count.table[,1]
      pep.count.table$count = pep.count.table$count+1
      
      # Log2 transform and impute missing values
      protein.matrix = log2(as.matrix(df.LFQ) + 1 )
      protein.matrix = impSeq(protein.matrix)
      protein.matrix = protein.matrix[,design.temp$File.Name]
      
      # Create design matrix for linear modeling
      design_model = model.matrix(~0+condition, data = design.temp)
      
      # Fit linear model and apply empirical Bayes moderation
      fit1 = lmFit(protein.matrix,design = design_model)
      cont <- makeContrasts(paste0("condition" , i.Contrasts1, "-condition" , i.Contrasts2), levels = design_model)
      fit2 = contrasts.fit(fit1,contrasts = cont)
      fit3 <- eBayes(fit2)
      
      # Assign peptide counts to model
      fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
      fit3$count[is.na(fit3$count)] <- 1
      
      # Check for valid peptide counts
      if(min(fit3$count) < 1){
        stop("fit3$count less than1")
      }
      
      # Apply DEqMS variance modeling
      fit4 = spectraCounteBayes(fit3)
      
      # Extract results
      DEqMS.results = outputResult(fit4,coef_col = 1)
      
      # Save results to file
      saveRDS(DEqMS.results, paste0("DIANN_results/", "DEqMS_" , i.exp , "_" , i.Contrasts, ".rds"))
      
      # Clean up memory
      remove(fit3)
      remove(fit4)
      remove(evidence_temp)
      gc()    
    }
  }
  
  # Clean up memory for next experiment
  remove(evidence)
  gc()
}
