# Load required libraries
library(limma)
library(stringr)
library(SummarizedExperiment)
library(LimROTS)
library(ROTS)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(samr)

# Load DIA-NN output for HEof without (n600) and with (w600) 600 proteins
exp.mar1 <- read.csv(paste0("../Spectronaut/", "HEof_n600" , "_DIA_Spectronaut_dlfq.tsv"), sep = "\t")
exp.mar1[exp.mar1 == 0] <- NA  # Replace 0s with NA for missing values
exp.mar2 <- read.csv(paste0("../Spectronaut/", "HEof_w600" , "_DIA_Spectronaut_dlfq.tsv"), sep = "\t")
exp.mar2[exp.mar2 == 0] <- NA

# Add suffixes to column names to identify batches
colnames(exp.mar1)[3:22] <- paste0(colnames(exp.mar1)[3:22] , "_N")
colnames(exp.mar2)[3:26] <- paste0(colnames(exp.mar2)[3:26] , "_W")

# Log2 transformation after imputation
exp.mar1[,c(3:22)] <- log2(exp.mar1[,c(3:22)]+ 1 )
exp.mar2[,c(3:26)] <- log2(exp.mar2[,c(3:26)]+ 1 )

# MDS plots for each batch
plotMDS(exp.mar1[,-c(1,2)])
plotMDS(exp.mar2[,-c(1,2)])

# Clean up protein and organism columns and set row names to protein ID
exp.mar1$Organism <- NULL
row.names(exp.mar1) <- exp.mar1$Protein
exp.mar1$Protein <- NULL

exp.mar2$Organism <- NULL
row.names(exp.mar2) <- exp.mar2$Protein
exp.mar2$Protein <- NULL

# Load and process design matrices for both batches
design1 <- read.csv(paste0("../Spectronaut/", "HEof_n600" , "_DIA_Spectronaut_design.tsv"), sep = "\t")
design2 <- read.csv(paste0("../Spectronaut/", "HEof_w600" , "_DIA_Spectronaut_design.tsv"), sep = "\t")

design1$sample_name <- paste0(design1$sample_name , "_N")
design2$sample_name  <- paste0(design2$sample_name , "_W")

design1$batch <- "N"
design2$batch <- "W"

# Combine design matrices
design <- rbind(design1, design2)

# Merge both expression matrices by protein ID
exp.mar <- merge(exp.mar1, exp.mar2, by = "row.names" , all = TRUE)
row.names(exp.mar) <- exp.mar$Row.names
exp.mar$Row.names <- NULL

# Filter columns to match those in the design matrix
exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]

# Ensure order matches between expression and design matrix
exp.mar <- exp.mar[,design$sample_name]

# Plot combined MDS
plotMDS(exp.mar)

# Generate pairwise condition contrasts
Contrasts <- combn(design$condition, 2, simplify = F)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]

# Reorder and filter again to match
exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]
exp.mar <- exp.mar[,design$sample_name]

# Main loop for analysis
i.exp <- "HEof"

for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    select.col.Contrasts1 <- nrow(design[design$condition %in% c(i.Contrasts1),])
    select.col.Contrasts2 <- nrow(design[design$condition %in% c(i.Contrasts2),])
    
    if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1){
        
        # Subset design and data for current contrast
        design.temp <- design[design$condition %in% c(i.Contrasts1, i.Contrasts2),]
        df.temp <- exp.mar[,colnames(exp.mar) %in% design.temp$sample_name]
        df.temp <- df.temp[,design.temp$sample_name]
        
        # Handle missing values
        if(any(is.na(df.temp))){
            df.temp <- impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        } else {
            df.temp[df.temp == 0 ] <- NA
            df.temp <- impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        }
        
        # Prepare metadata
        sample_info <- data.frame(
            sample_id = colnames(df.temp),
            group = design.temp$condition,
            batches = design.temp$batch)
        
        rownames(sample_info) <- row.names(sample_info$sample_id)
        sample_info$group <- as.factor(as.numeric(as.factor(sample_info$group)))
        sample_info$batches <- as.factor(sample_info$batches)
        
        feature_info <- data.frame(protein_id = row.names(df.temp))
        rownames(feature_info) <- feature_info$protein_id
        
        # Create SummarizedExperiment
        se <- SummarizedExperiment(
            assays = list(protein.exp = df.temp),
            colData = sample_info,
            rowData = feature_info 
        )
        
        # Run LimROTS
        meta.info <- c("group" , "batches")
        niter <- 500
        K <- floor(nrow(feature_info) / 2)
        group.name <- "group"
        formula.str <- "~0+group+group:batches"
        
        se <- LimROTS(
            x = se,
            niter = niter, K = K, meta.info = meta.info,
            BPPARAM  = SnowParam(5, progressbar = T), group.name = group.name,
            formula.str = formula.str, trend = TRUE,
            robust = TRUE, permutating.group = FALSE
        )
        
        saveRDS(se, paste0("Spectronaut_results/", "LimROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(se); gc()
        
        # Run ROTS
        rots_results <- ROTS(data = df.temp , groups = as.numeric(sample_info$group), B = niter, K = K, 
                             seed = 1597, progress = TRUE, verbose = TRUE)
        saveRDS(rots_results, paste0("Spectronaut_results/", "ROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(rots_results); gc()
        
        # Run Limma
        sample_info$group <- as.factor(design.temp$condition) 
        design.model = model.matrix(~0+group+group:batches , data = sample_info) 
        colnames(design.model) <- make.names(colnames(design.model))
        
        fit1 = lmFit(df.temp, design = design.model)
        cont <- makeContrasts(contrasts =  paste0("group", i.Contrasts1 , "-" , "group" , i.Contrasts2), levels = design.model)
        fit2 = contrasts.fit(fit1,contrasts = cont)
        fit3 <- eBayes(fit2, trend = T, robust = T)
        limma.results = topTable(fit3, adjust="BH", sort.by = 'logFC', n=Inf)
        saveRDS(limma.results, paste0("Spectronaut_results/", "Limma_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(limma.results); gc()
        
        # Run SAM
        data <- list(x=as.matrix(df.temp), y=sample_info$group,
                     geneid=feature_info$protein_id,
                     genenames=feature_info$protein_id, logged2=TRUE)
        
        samr.obj <- samr(data, resp.type="Two class unpaired", nperms=100)
        logFC <- log2(samr.obj$foldchange)
        pvalue <- data.frame(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
        adj.pvalue <- p.adjust(pvalue[,1], method = "BH")
        SAM.results <- cbind(logFC, pvalue, adj.pvalue)
        saveRDS(SAM.results, paste0("Spectronaut_results/", "SAM_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(SAM.results); gc()
        
        # Run t-test
        i.test.result <- data.frame()
        for(i.test in seq_len(nrow(df.temp))){
            i.test.df <- df.temp[i.test,]
            pvalue <- t.test(formula = as.numeric(i.test.df[1,]) ~ sample_info$group)$p.value
            fc.calc <- mean(as.numeric(i.test.df[,1:select.col.Contrasts1])) -
                       mean(as.numeric(i.test.df[,select.col.Contrasts1+1:select.col.Contrasts2]))
            
            i.test.result <- rbind(i.test.result, data.frame(row.names = row.names(i.test.df), 
                                                             pvalue = pvalue , 
                                                             fc.calc = fc.calc))
        }
        i.test.result$adj.pvalue <- p.adjust(i.test.result$pvalue, method = "BH")
        saveRDS(i.test.result, paste0("Spectronaut_results/", "ttest_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(i.test.result); gc()
        
        # run ANOVA 

        design.model = model.matrix(~0+group+group:batches , data = sample_info) 
        
        fit1 = lmFit(df.temp,
                     design = design.model)
        ord.t = fit1$coefficients[, 1]/fit1$sigma/fit1$stdev.unscaled[, 1]
        ord.p = 2*pt(abs(ord.t), fit1$df.residual, lower.tail = FALSE)
        ord.q = p.adjust(ord.p,method = "BH")
        log2FC = fit1$coefficients[,2]-fit1$coefficients[,1]
        anova.results = data.frame(row.names(df.temp),
                                   logFC=log2FC,
                                   t=ord.t, 
                                   P.Value=ord.p, 
                                   adj.P.Val = ord.q)
        saveRDS(anova.results, paste0("Spectronaut_results/", "ANOVA_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(anova.results)
        gc()

        # run DEP

        experimental_design <- design.temp[,c(3,4,5,6)]
        colnames(experimental_design)[1] <- "label"
        df.temp$name <- row.names(df.temp)
        df.temp$ID <- row.names(df.temp)
        experimental_design$replicate <- seq(1, nrow(experimental_design))
        
        
        
        data_se <- make_se(df.temp, seq(1, ncol(df.temp)-2), experimental_design)
        
        data_se$label <- make.names(data_se$label)
        data_se$condition <- make.names(data_se$condition)
        data_se$batch <- make.names(data_se$batch)
        data_se$batch <- as.factor(data_se$batch)
        
        
        data_diff_manual <- test_diff(data_se, control = i.Contrasts1 , test = i.Contrasts2,
                                      design_formula = formula("~ 0 + condition+batch"))
        
        data_diff_manual.df <- elementMetadata(data_diff_manual)
        
        data_diff_manual.df <- data.frame(data_diff_manual.df)
        
        data_diff_manual.df <- data_diff_manual.df[,c(5,6,7)]
        colnames(data_diff_manual.df) <- str_remove(colnames(data_diff_manual.df) , paste0(i.Contrasts2 , "_vs_" , i.Contrasts1) )
        colnames(data_diff_manual.df) <- str_remove(colnames(data_diff_manual.df) , fixed("_") )
        
        
        
        saveRDS(data_diff_manual.df, paste0("Spectronaut_results/", "DEP_" , i.exp , "_" , i.Contrasts, ".rds"))
        
        
    }
}
