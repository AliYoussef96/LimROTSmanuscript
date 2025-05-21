# Load required libraries
library(stringr)     # For string manipulation
library(MSstats)     # For mass spectrometry data analysis
library(DEqMS)       # For protein-level differential expression adjusted by peptide count
library(dplyr)       # For data manipulation
library(rrcovNA)     # For robust estimation with missing values

# Set experiment identifier
i.exp <- "HEof"

# Load narrow-window DIA data and its annotation
exp.mar1 <- read.csv("../HEof_n600_DIA/Spectronaut/20231013_143752_Narrow_Report.tsv", sep = "\t")
annot.mar1 <- read.csv("../Spectronaut/HEof_n600_DIA_spt_design_msstats.tsv", sep = "\t")

# Add suffix to distinguish conditions from narrow window
annot.mar1$Condition <- paste0(annot.mar1$Condition , "_Narrow")
exp.mar1$R.Condition <- paste0(exp.mar1$R.Condition , "_Narrow")

# Load wide-window DIA data and its annotation
exp.mar2 <- read.csv("../HEof_w600_DIA/Spectronaut/20231013_160348_25fmol_Report.tsv", sep = "\t")
annot.mar2 <- read.csv("../Spectronaut/HEof_w600_DIA_spt_design_msstats.tsv", sep = "\t")

# Add suffix to distinguish conditions from wide window
annot.mar2$Condition <- paste0(annot.mar2$Condition , "_Wide")
exp.mar2$R.Condition <- paste0(exp.mar2$R.Condition , "_Wide")

# Combine narrow and wide experiments and annotations
exp <- rbind(exp.mar1, exp.mar2)
annot <- rbind(annot.mar1, annot.mar2)

############################

# Load protein-level intensities (DLFQ) for narrow and wide datasets
exp.mar1 <- read.csv(paste0("../Spectronaut/", "HEof_n600" , "_DIA_spt_dlfq.tsv"), sep = "\t")
exp.mar1[exp.mar1 == 0] <- NA   # Convert 0 to NA

exp.mar2 <- read.csv(paste0("../Spectronaut/", "HEof_w600" , "_DIA_spt_dlfq.tsv"), sep = "\t")
exp.mar2[exp.mar2 == 0] <- NA   # Convert 0 to NA

# Rename columns to indicate batch
colnames(exp.mar1)[3:25] <- paste0(colnames(exp.mar1)[3:25] , "_N")
colnames(exp.mar2)[3:26] <- paste0(colnames(exp.mar2)[3:26] , "_W")

# Log2-transform intensity values (+1 to avoid log(0))
exp.mar1[,c(3:25)] <- log2(exp.mar1[,c(3:25)]+ 1 )
exp.mar2[,c(3:26)] <- log2(exp.mar2[,c(3:26)]+ 1 )

# Prepare for merging: set row names and remove unnecessary columns
exp.mar1$Organism <- NULL
row.names(exp.mar1) <- exp.mar1$Protein
exp.mar1$Protein <- NULL

exp.mar2$Organism <- NULL
row.names(exp.mar2) <- exp.mar2$Protein
exp.mar2$Protein <- NULL

# Load design files for both batches
design1 <- read.csv(paste0("../Spectronaut/", "HEof_n600" , "_DIA_spt_design.tsv"), sep = "\t")
design2 <- read.csv(paste0("../Spectronaut/", "HEof_w600" , "_DIA_spt_design.tsv"), sep = "\t")

# Adjust sample names and assign batch labels
design1$sample_name <- paste0(design1$sample_name , "_N")
design2$sample_name  <- paste0(design2$sample_name , "_W")
design1$batch <- "N"
design2$batch <- "W"

# Merge design matrices
design <- rbind(design1, design2)

# Merge protein expression data by row names
exp.mar <- merge(exp.mar1, exp.mar2, by = "row.names" , all = TRUE)
row.names(exp.mar) <- exp.mar$Row.names

# Subset expression data to match samples in the design matrix
exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]
exp.mar <- exp.mar[,design$sample_name]

###########
# Define pairwise comparisons for differential testing
Contrasts <- combn(c("A", "B" , "C" , "D" , "E" , "F", "G" , "H"), 2, simplify = F)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))  # Sort and collapse
Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]  # Remove same-group comparisons

# Set random seed for reproducibility
set.seed(123456, sample.kind = "default" , kind = "default")

# Loop through each contrast
for(i.Contrasts in Contrasts){
  print(i.Contrasts)
  
  # Extract condition labels for contrast
  i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
  i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
  
  # Create label names for narrow and wide versions
  i.Contrasts1.n <- paste0(i.Contrasts1 , "_Narrow")
  i.Contrasts1.w <- paste0(i.Contrasts1 , "_Wide")
  i.Contrasts2.n <- paste0(i.Contrasts2 , "_Narrow")
  i.Contrasts2.w <- paste0(i.Contrasts2 , "_Wide")
  
  # Subset design matrix for current comparison
  design.temp <- design[design$condition %in% c(i.Contrasts1,i.Contrasts2 ),]
  df.prot <- exp.mar
  design.temp <- design.temp[order(design.temp$condition),]
  df.prot <- df.prot[,design.temp$sample_name]
  
  # Impute missing values using sequential imputation
  df.prot.filter <- impSeq(df.prot)
  
  # Filter peptide quantification table to only those used in this contrast
  exp_ <- exp[exp$R.FileName %in% design.temp$file,]
  
  # Count unique peptides per protein for DEqMS modeling
  pep.count.table_ <- exp_[!duplicated(exp[,c(7,13)]),]
  pep.count.table_ <- pep.count.table_ %>% group_by(PG.ProteinAccessions) %>% summarise(n = n())
  pep.count.table_ <- data.frame(pep.count.table_)
  pep.count.table_ <- pep.count.table_[!is.na(pep.count.table_$n),]
  
  # Match protein IDs in expression data to those in the count table
  pep.count.table <- data.frame()
  for(i in unique(row.names(df.prot.filter))){
    i_ <- ifelse(grepl("ECOLI" , i),
                 stringr::str_extract_all(i, "(?<=sp\\|)[^|]+")[[1]] |> paste(collapse = ";"),
                 stringr::str_extract_all(i, "(?<=\\|)[^|]+")[[1]][seq(2, by = 2, length.out = length(stringr::str_extract_all(i, "(?<=\\|)[^|]+")[[1]]) / 2)] |> paste(collapse = ";"))
    temp <- pep.count.table_[pep.count.table_$PG.ProteinAccessions == i_,]$n
    pep.count.table <- rbind(pep.count.table , data.frame(protein = i, count = temp))
  }
  
  pep.count.table <- pep.count.table[!is.na(pep.count.table$count),]
  pep.count.table$count = pep.count.table$count + 1  # Add pseudocount to avoid 0
  
  # Fit linear model with condition and condition:batch interaction
  design_model = model.matrix(~0+condition+condition:batch, data = design.temp)
  print(design_model)
  
  # Fit model using limma
  df.prot.filter <- df.prot.filter[,design.temp$sample_name]
  colnames(design_model) <- make.names(colnames(design_model))
  fit1 = lmFit(df.prot.filter, design = design_model)
  
  # Create contrast for the specific condition comparison
  cont <- makeContrasts(paste0("condition" , i.Contrasts1, "-condition" , i.Contrasts2), levels = design_model)
  fit2 = contrasts.fit(fit1, contrasts = cont)
  fit3 <- eBayes(fit2)
  
  # Add peptide count to model fit
  row.names(pep.count.table) <- pep.count.table$protein
  fit3$count = pep.count.table[rownames(fit3$coefficients), "count"]
  fit3$count <- ifelse(is.na(fit3$count), 1 , fit3$count)  # Replace missing with 1
  
  # Stop if any count < 1 (shouldn’t happen due to pseudocount)
  if(min(fit3$count) < 1){
    stop("fit3$count less than 1")
  }
  
  # Adjust variance using peptide count
  fit4 = spectraCounteBayes(fit3)
  
  # Get DEqMS results
  DEqMS.results = outputResult(fit4, coef_col = 1)
  
  # Save results as RDS file
  saveRDS(DEqMS.results, paste0("results_case2_DIAN/", "DEqMS_" , i.exp , "_" , i.Contrasts, ".rds"))
}
