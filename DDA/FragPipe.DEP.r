library(stringr)
library(SummarizedExperiment)
library(LimROTS)
library(ROTS)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(limma)
library(samr)
library(multiUS)

# List all files in the directory and extract experiment names
all.exp  <- list.files("../FragPipe/")
all.exp <- str_split_fixed(all.exp, "_LFQ_FragPipe_" , 2)[,1]
all.exp <- unique(all.exp)

# Set random seed
set.seed(1597, sample.kind = "default" , kind = "default")

# Loop through each experiment
for(i.exp in all.exp){
    print(i.exp)
    
    # Read intensity and design files
    exp.mar <- read.csv(paste0("../FragPipe/", i.exp , "_LFQ_FragPipe_dlfq_pro_intensity.tsv"), sep = "\t")
    design <- read.csv(paste0("../FragPipe/", i.exp , "_LFQ_FragPipe_design.tsv"), sep = "\t")
    
    # Generate pairwise contrasts of conditions
    Contrasts <- combn(design$condition, 2, simplify = F)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
    
    # Format intensity data
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    # Subset matching columns and samples
    exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
    design <- design[design$sample_name %in% colnames(exp.mar),]
    exp.mar <- exp.mar[,design$sample_name]
    
    # Loop through each valid contrast
    for(i.Contrasts in Contrasts){
        print(i.Contrasts)
        
        # Extract two conditions from contrast
        i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
        i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
        
        # Count replicates per condition
        select.col.Contrasts1 <- nrow(design[design$condition %in% c(i.Contrasts1),])
        select.col.Contrasts2 <- nrow(design[design$condition %in% c(i.Contrasts2),])
        
        # Proceed only if both conditions have more than 1 replicate
        if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1){
            
            # Subset design and data for current contrast
            design.temp <- design[design$condition %in% c(i.Contrasts1, i.Contrasts2),]
            df.temp <- exp.mar[,colnames(exp.mar) %in% design.temp$sample_name]
            df.temp <- df.temp[,design.temp$sample_name]
            
            # Impute missing values
    
            df.temp[df.temp == 0 ] <- NA
            df.temp <- log2(df.temp+1)
            df.temp <- seqKNNimp(df.temp)
            df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        
            
            # Prepare experimental design
            experimental_design <- design.temp[,c(3,4,5)]
            colnames(experimental_design)[1] <- "label"
            df.temp$name <- row.names(df.temp)
            df.temp$ID <- row.names(df.temp)
            experimental_design$replicate <- seq(1, nrow(experimental_design))
            
            # Create SummarizedExperiment object
            data_se <- make_se(df.temp, seq(1, ncol(df.temp)-2), experimental_design)
            
            # Clean condition and label names
            data_se$label <- make.names(data_se$label)
            data_se$condition <- make.names(data_se$condition)
            
            # Run differential analysis
            data_diff_manual <- test_diff(data_se, control = i.Contrasts1 , test = i.Contrasts2,
                                          design_formula = formula("~ 0 + condition"))
            
            # Extract and format results
            data_diff_manual.df <- elementMetadata(data_diff_manual)
            data_diff_manual.df <- data.frame(data_diff_manual.df)
            data_diff_manual.df <- data_diff_manual.df[,c(5,6,7)]
            colnames(data_diff_manual.df) <- str_remove(colnames(data_diff_manual.df) , paste0(i.Contrasts2 , "_vs_" , i.Contrasts1) )
            colnames(data_diff_manual.df) <- str_remove(colnames(data_diff_manual.df) , fixed("_") )
            
            # Save results
            saveRDS(data_diff_manual.df, paste0("FragPipe_results/", "DEP_" , i.exp , "_" , i.Contrasts, ".rds"))
        }
    }
}
