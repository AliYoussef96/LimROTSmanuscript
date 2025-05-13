library(limma)
library(stringr)
library(SummarizedExperiment)
library(LimROTS)
library(ROTS)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(limma)
library(samr)


exp.mar1 <- read.csv(paste0("../DIANN/", "HEof_n600" , "_DIA_DIANN_dlfq.tsv"), sep = "\t")
exp.mar1[exp.mar1 == 0] <- NA
exp.mar2 <- read.csv(paste0("../DIANN/", "HEof_w600" , "_DIA_DIANN_dlfq.tsv"), sep = "\t")
exp.mar2[exp.mar2 == 0] <- NA

colnames(exp.mar1)[3:22] <- paste0(colnames(exp.mar1)[3:22] , "_N")
colnames(exp.mar2)[3:26] <- paste0(colnames(exp.mar2)[3:26] , "_W")


exp.mar1[,c(3:22)] <- log2(exp.mar1[,c(3:22)]+ 1 )
exp.mar2[,c(3:26)] <- log2(exp.mar2[,c(3:26)]+ 1 )

plotMDS(exp.mar1[,-c(1,2)])
plotMDS(exp.mar2[,-c(1,2)])

exp.mar1$Organism <- NULL
row.names(exp.mar1) <- exp.mar1$Protein
exp.mar1$Protein <- NULL

exp.mar2$Organism <- NULL
row.names(exp.mar2) <- exp.mar2$Protein
exp.mar2$Protein <- NULL

design1 <- read.csv(paste0("../DIANN/", "HEof_n600" , "_DIA_DIANN_design.tsv"), sep = "\t")
design2 <- read.csv(paste0("../DIANN/", "HEof_w600" , "_DIA_DIANN_design.tsv"), sep = "\t")

design1$sample_name <- paste0(design1$sample_name , "_N")
design2$sample_name  <- paste0(design2$sample_name , "_W")

design1$batch <- "N"
design2$batch <- "W"

design <- rbind(design1, design2)

exp.mar <- merge(exp.mar1, exp.mar2, by = "row.names" , all = TRUE)

row.names(exp.mar) <- exp.mar$Row.names
exp.mar$Row.names <- NULL

exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]

exp.mar <- exp.mar[,design$sample_name]

plotMDS(exp.mar)

###########

Contrasts <- combn(design$condition, 2, simplify = F)
Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]

#exp.mar$Organism <- NULL
#design$sample_name <- str_remove(design$sample_name, fixed("_"))
#colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))

exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
design <- design[design$sample_name %in% colnames(exp.mar),]

exp.mar <- exp.mar[,design$sample_name]

i.exp <- "HEof"

for(i.Contrasts in Contrasts){
    print(i.Contrasts)
    
    i.Contrasts1 <- str_split(i.Contrasts, "")[[1]][1]
    i.Contrasts2 <- str_split(i.Contrasts, "")[[1]][2]
    
    
    select.col.Contrasts1 <- nrow(design[design$condition %in% c(i.Contrasts1),])
    select.col.Contrasts2 <- nrow(design[design$condition %in% c(i.Contrasts2),])
    
    if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1){
        
        design.temp <- design[design$condition %in% c(i.Contrasts1, i.Contrasts2),]
        
        df.temp <- exp.mar[,colnames(exp.mar) %in% design.temp$sample_name]
        
        df.temp <- df.temp[,design.temp$sample_name]
        
        if(any(is.na(df.temp))){
            #df.temp <- log2(df.temp+1)
            df.temp <-  impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        }else{
            #df.temp[df.temp == 0 ] <- NA
            df.temp <- log2(df.temp+1)
            df.temp <-  impute.MinDet(df.temp)
            df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        }
        
        add.effect <- design.temp[design.temp$batch == "W" & design.temp$condition == design.temp$condition[1],]
        ecoli <- sample(which(grepl("ECOLI", row.names(df.temp))), 500)
        
        df.temp[ecoli,which(colnames(df.temp) %in% add.effect$sample_name)] <- df.temp[ecoli,which(colnames(df.temp) %in% add.effect$sample_name)]  + runif(500 , min = 5, max = 20)
        
        row.names(df.temp)[ecoli] <- paste0(row.names(df.temp)[ecoli] , "_" , "added")
        
        sample_info <- data.frame(
            sample_id = colnames(df.temp),
            group = design.temp$condition,
            batches = design.temp$batch)
        
        rownames(sample_info) <- row.names(sample_info$sample_id)
        sample_info$group <- as.factor( as.numeric(as.factor(sample_info$group)) )
        sample_info$batches <- as.factor(sample_info$batches)
        
        feature_info <- data.frame(
            protein_id = row.names(df.temp)
        )
        rownames(feature_info) <- feature_info$protein_id
        
        
        # Create the SummarizedExperiment object
        se <- SummarizedExperiment(
            assays = list(protein.exp = df.temp),
            colData = sample_info,
            rowData = feature_info 
        )
        
        # Set metadata and formula for LimROTS analysis
        meta.info <- c("group" , "batches")
        niter <- 500 # Number of bootstrap samples
        K <- nrow(feature_info) / 2 # Set the value for K based on the data size
        K <- floor(K)
        group.name <- "group"
        formula.str <- "~0+group+group:batches" # Formula for group comparison
        
        # Run LimROTS analysis with trend and robust settings enabled
        se <- LimROTS(
            x = se,
            niter = niter, K = K, meta.info = meta.info,
            BPPARAM  = SnowParam(5, progressbar = T), group.name = group.name,
            formula.str = formula.str, trend = TRUE,
            robust = TRUE, permutating.group = FALSE
        )
        
        saveRDS(se, paste0("DIANN_results/", "LimROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
        
        remove(se)
        gc()
        
        rots_results <- ROTS(data = df.temp , groups = as.numeric(sample_info$group), B = niter, K = K, 
                             seed = 1597, progress = TRUE, verbose = TRUE)
        
        saveRDS(rots_results, paste0("DIANN_results/", "ROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
        
        remove(rots_results)
        gc()
        
        sample_info$group <- as.factor(design.temp$condition) 
        
        design.model = model.matrix(~0+group+group:batches , data = sample_info) 
        colnames(design.model) <- make.names(colnames(design.model))
        
        fit1 = lmFit(df.temp, design = design.model)
        cont <- makeContrasts(contrasts =  paste0( "group", i.Contrasts1 , "-" , "group" , i.Contrasts2), levels = design.model)
        fit2 = contrasts.fit(fit1,contrasts = cont)
        fit3 <- eBayes(fit2, trend = T, robust = T)
        limma.results = topTable(fit3, adjust="BH", sort.by = 'logFC', n=Inf)
        
        saveRDS(limma.results, paste0("DIANN_results/", "Limma_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(limma.results)
        gc()
        
        
        
        # run SAM
        
        data<-list(x=as.matrix(df.temp), y=sample_info$group
                   , geneid=feature_info$protein_id,
                   genenames=feature_info$protein_id, logged2=TRUE)
        
        samr.obj <- samr( data, resp.type="Two class unpaired", nperms=100)
        logFC<-log2(samr.obj$foldchange)
        
        pvalue <- data.frame(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
        
        adj.pvalue <- p.adjust(pvalue$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar., method = "BH")
        
        SAM.results<-cbind(logFC, pvalue, adj.pvalue)
        
        saveRDS(SAM.results, paste0("DIANN_results/", "SAM_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(SAM.results)
        gc()
        
        
        # run t-test
        i.test.result <- data.frame()
        for(i.test in seq_len(nrow(df.temp))){
            i.test.df <- df.temp[i.test,]
            pvalue <- t.test(formula = as.numeric(i.test.df[1,])~sample_info$group)$p.value
            fc.calc <- mean( as.numeric( i.test.df[,1:select.col.Contrasts1] ) ) - mean( as.numeric(i.test.df[,select.col.Contrasts1+1:select.col.Contrasts2]) )
            
            i.test.result <- rbind(i.test.result, data.frame(row.names = row.names(i.test.df), 
                                                             pvalue = pvalue , 
                                                             fc.calc = fc.calc))
        }
        
        i.test.result$adj.pvalue <- p.adjust(i.test.result$pvalue, method = "BH")
        
        
        saveRDS(i.test.result, paste0("DIANN_results/", "ttest_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(i.test.result)
        gc()
        
        # run ANOVA
        
        
        ## codes from:
        #https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#compare-p-values-from-deqms-to-ordinary-t-test-anova-and-limma
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
        saveRDS(anova.results, paste0("DIANN_results/", "ANOVA_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(anova.results)
        gc()
        
        
        }
    
    }




