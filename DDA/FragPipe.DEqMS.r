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
library(multiUS)

# List all experiment folders except specific ones
all.exp  <- list.files("../../")
all.exp <- all.exp[!all.exp %in% c("FragPipe" , "case study 1", "Maxquant")]

# Set reproducible seed
set.seed(1597, sample.kind = "default" , kind = "default")

# Loop over all experiment folders
for(i.exp in all.exp){
    print(i.exp)
    
    # Read protein quantification data
    exp.mar <- read.csv(paste0("../../", i.exp, "/FragPipe/TOP0/combined_protein.tsv"  ), sep = "\t")
    
    # Read design file
    design <- read.csv(paste0("../../FragPipe/", i.exp , "_FragPipe_design.tsv"), sep = "\t")
    
    # Read annotation for MSstats (not used below)
    annot <- read.csv(paste0("../../FragPipe/", i.exp , "_FragPipe_design_msstats.tsv"), sep = "\t")
    
    # Define unique pairwise contrasts
    Contrasts <- combn(design$condition, 2, simplify = F)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
    
    # Format data
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Loop over all valid contrasts
    for(i.Contrasts in Contrasts){
        print(i.Contrasts)
        
        # Get contrast groups
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        # Subset design to relevant samples
        design.temp <- design[design$condition %in% c(i.Contrasts1,i.Contrasts2 ),]
        design.temp$LFQ.samples <- paste0(design.temp$sample_name, ".MaxLFQ.Intensity")
        
        # Read protein data again
        df.prot <- read.csv(paste0("../../", i.exp, "/FragPipe/TOP0/combined_protein.tsv"  ), sep = "\t")
        
        # Extract LFQ columns and replace 0s with NA
        df.LFQ = df.prot[,which(grepl("MaxLFQ", colnames(df.prot)))]
        df.LFQ[df.LFQ==0] <- NA
        
        # Format LFQ matrix
        rownames(df.LFQ) = df.prot$Protein
        colnames(df.LFQ) <- str_remove(colnames(df.LFQ)  , fixed("_"))
        df.LFQ <- df.LFQ[,colnames(df.LFQ) %in% design.temp$LFQ.samples]
        design.temp <- design.temp[order(design.temp$condition),]
        
        # Ensure all design samples are in LFQ matrix
        if(! any( !design.temp$LFQ.samples %in% colnames(df.LFQ) ))
        {
            df.LFQ <- df.LFQ[,design.temp$LFQ.samples]
            
            # Create peptide count table and add pseudocount
            pep.count.table = data.frame(count = df.prot$Combined.Total.Peptides,
                                         row.names = df.prot$Protein)
            pep.count.table$count = pep.count.table$count+1
            
            # Log2 transform LFQ matrix with pseudocount
            protein.matrix = log2(as.matrix(df.LFQ) + 1 )
            
            # Impute missing values using KNN
            protein.matrix = seqKNNimp(protein.matrix)
            
            # Build design matrix for limma model
            design_model = model.matrix(~0+condition, data = design.temp)
            print(design_model)
            
            # Fit limma model and apply contrast
            fit1 = lmFit(protein.matrix,design = design_model)
            cont <- makeContrasts(paste0("condition" , i.Contrasts1, "-condition" , i.Contrasts2), levels = design_model)
            fit2 = contrasts.fit(fit1,contrasts = cont)
            fit3 <- eBayes(fit2)
            
            # Add peptide count to fit object
            fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
            
            # Check for invalid counts
            if(min(fit3$count) < 1){
                stop("fit3$count less than1")
            }
            
            # Apply DEqMS with peptide count adjustment
            fit4 = spectraCounteBayes(fit3)
            
            # Format DEqMS output
            DEqMS.results = outputResult(fit4,coef_col = 1)
            rownames(df.prot) = df.prot$Majority.protein.IDs
            DEqMS.results$Gene.name = df.prot[DEqMS.results$gene,]$Gene.names
            
            # Save results to RDS
            saveRDS(DEqMS.results, paste0("FragPipe_results/", "DEqMS_" , i.exp , "_" , i.Contrasts, ".rds"))
            
            # Cleanup
            remove(fit3)
            remove(fit4)
            gc()
        }
    }
}
