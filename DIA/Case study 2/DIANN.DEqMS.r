# Load required libraries
library(stringr)
library(MSstats)
library(DEqMS)
library(dplyr)
library(rrcovNA)

# Set experiment name
i.exp <- "HEof"

# Read narrow window DIA experiment data and annotation
exp.mar1 <- read.csv("../HEof_n600_DIA/DIANN/report.tsv", sep = "\t")
annot.mar1 <- read.csv("../DIANN/HEof_n600_DIA_DIANN_design_msstats.tsv", sep = "\t")
annot.mar1$Condition <- paste0(annot.mar1$Condition , "_Narrow")  # Append "_Narrow" to condition
exp.mar1$R.Condition <- paste0(exp.mar1$R.Condition , "_Narrow")  # Append "_Narrow" to reported condition

# Read wide window DIA experiment data and annotation
exp.mar2 <- read.csv("../HEof_w600_DIA/DIANN/report.tsv", sep = "\t")
annot.mar2 <- read.csv("../DIANN/HEof_w600_DIA_DIANN_design_msstats.tsv", sep = "\t")
annot.mar2$Condition <- paste0(annot.mar2$Condition , "_Wide")  # Append "_Wide" to condition
exp.mar2$R.Condition <- paste0(exp.mar2$R.Condition , "_Wide")  # Append "_Wide" to reported condition

# Combine experiments and annotations
exp <- rbind(exp.mar1, exp.mar2)
annot <- rbind(annot.mar1, annot.mar2)

############################

# Read DIA LFQ quantification files and replace 0s with NA
exp.mar1 <- read.csv(paste0("../DIANN/", "HEof_n600" , "_DIA_DIANN_dlfq.tsv"), sep = "\t")
exp.mar1[exp.mar1 == 0] <- NA
exp.mar2 <- read.csv(paste0("../DIANN/", "HEof_w600" , "_DIA_DIANN_dlfq.tsv"), sep = "\t")
exp.mar2[exp.mar2 == 0] <- NA

# Rename columns to reflect narrow (N) and wide (W) batches
colnames(exp.mar1)[3:22] <- paste0(colnames(exp.mar1)[3:22] , "_N")
colnames(exp.mar2)[3:26] <- paste0(colnames(exp.mar2)[3:26] , "_W")

# Log2 transformation of intensities (+1 to avoid log(0))
exp.mar1[,c(3:22)] <- log2(exp.mar1[,c(3:22)]+ 1 )
exp.mar2[,c(3:26)] <- log2(exp.mar2[,c(3:26)]+ 1 )

# Remove metadata columns and set row names to protein identifiers
exp.mar1$Organism <- NULL
row.names(exp.mar1) <- exp.mar1$Protein
exp.mar1$Protein <- NULL

exp.mar2$Organism <- NULL
row.names(exp.mar2) <- exp.mar2$Protein
exp.mar2$Protein <- NULL

# Read sample design files and label batch
design1 <- read.csv(paste0("../DIANN/", "HEof_n600" , "_DIA_DIANN_design.tsv"), sep = "\t")
design2 <- read.csv(paste0("../DIANN/", "HEof_w600" , "_DIA_DIANN_design.tsv"), sep = "\t")
design1$sample_name <- paste0(design1$sample_name , "_N")
design2$sample_name  <- paste0(design2$sample_name , "_W")
design1$batch <- "N"
design2$batch <- "W"
design <- rbind(design1, design2)

# Merge narrow and wide LFQ data
exp.mar <- merge(exp.mar1, exp.mar2, by = "row.names" , all = TRUE)
row.names(exp.mar) <- exp.mar$Row.names

# Keep only samples that are in the design table
exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]
exp.mar <- exp.mar[,design$sample_name]

# Generate pairwise contrasts of all conditions (excluding self-comparisons)
Contrasts <- combn(c("A", "B" , "C" , "D" , "E" , "F", "G" , "H"), 2, simplify = F)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]

# Set seed for reproducibility
set.seed(123456, sample.kind = "default" , kind = "default")

# Loop over all pairwise condition contrasts
for(i.Contrasts in Contrasts){
  print(i.Contrasts)
  
  # Parse conditions from contrast string
  i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
  i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
  
  # Label with batch suffixes
  i.Contrasts1.n <- paste0(i.Contrasts1 , "_Narrow")
  i.Contrasts1.w <- paste0(i.Contrasts1 , "_Wide")
  i.Contrasts2.n <- paste0(i.Contrasts2 , "_Narrow")
  i.Contrasts2.w <- paste0(i.Contrasts2 , "_Wide")
  
  # Subset design for the current contrast
  design.temp <- design[design$condition %in% c(i.Contrasts1,i.Contrasts2 ),]
  df.prot <- exp.mar
  design.temp <- design.temp[order(design.temp$condition),]
  df.prot <- df.prot[,design.temp$sample_name]
  
  # Impute missing values (sequential method)
  df.prot.filter <- impSeq(df.prot)
  
  # Filter corresponding raw file data for peptide count analysis
  exp_ <- exp[exp$File.Name %in% design.temp$file,]
  
  # Count unique peptides per protein
  pep.count.table_ <- exp_[!duplicated(exp_[,c(3,15)]),]
  pep.count.table_$ids <- ifelse(!grepl("HUMAN", pep.count.table_$Protein.Ids), paste0(pep.count.table_$Protein.Group, "|" , pep.count.table_$Protein.Names), pep.count.table_$Protein.Group)
  pep.count.table_ <- pep.count.table_ %>% group_by(ids) %>% summarise(n = n())
  pep.count.table <- data.frame()
  pep.count.table_ <- data.frame(pep.count.table_)
  
  # Map peptide counts to protein table
  for(i in unique(row.names(df.prot.filter))){
    i_ <- str_remove(i , fixed("sp|"))
    temp <- pep.count.table_[pep.count.table_$ids == i_,]$n
    if(!any(pep.count.table_$ids == i_)){
      temp <- 0
    }
    pep.count.table <- rbind(pep.count.table , data.frame(protein = i, count = temp))
  }

  # Add pseudocount to avoid zero values
  pep.count.table$count = pep.count.table$count+1
  pep.count.table <- pep.count.table %>% group_by(protein) %>% summarise(count = mean(count))
  pep.count.table <- data.frame(pep.count.table)

  # Create linear model design matrix (no intercept)
  design_model = model.matrix(~0+condition+condition:batch, data = design.temp)
  print(design_model)
  
  # Ensure sample order matches
  df.prot.filter <- df.prot.filter[,design.temp$sample_name]
  colnames(design_model) <- make.names(colnames(design_model))
  
  # Fit linear model with limma
  fit1 = lmFit(df.prot.filter,design = design_model)
  cont <- makeContrasts(paste0("condition" , i.Contrasts1, "-condition" , i.Contrasts2), levels = design_model)
  fit2 = contrasts.fit(fit1,contrasts = cont)
  fit3 <- eBayes(fit2)
  
  # Add peptide counts for DEqMS
  row.names(pep.count.table) <- pep.count.table$protein
  fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]

  # Ensure all proteins have at least one peptide
  if(min(fit3$count) < 1){
    stop("fit3$count less than1")
  }
  
  # Perform DEqMS analysis
  fit4 = spectraCounteBayes(fit3)
  DEqMS.results = outputResult(fit4,coef_col = 1)
  
  # Save results
  saveRDS(DEqMS.results, paste0("results_case2_DIAN/", "DEqMS_" , i.exp , "_" , i.Contrasts, ".rds"))
}
