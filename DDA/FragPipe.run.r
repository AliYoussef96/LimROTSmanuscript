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

all.exp  <- list.files("FragPipe/")
all.exp <- str_split_fixed(all.exp, "_LFQ_FragPipe_" , 2)[,1]
all.exp <- unique(all.exp)

set.seed(1597, sample.kind = "default" , kind = "default")


#### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### 
#### LimROTS, ROTS, Limma, SAM, ANOVA and t-test
#### #### #### #### #### #### #### #### #### 
#### #### #### #### #### #### #### #### #### 


for(i.exp in all.exp[10:12]){
    print(i.exp)
    
    exp.mar <- read.csv(paste0("FragPipe/", i.exp , "_LFQ_FragPipe_dlfq_pro_intensity.tsv"), sep = "\t")
    design <- read.csv(paste0("FragPipe/", i.exp , "_LFQ_FragPipe_design.tsv"), sep = "\t")
    Contrasts <- combn(design$condition, 2, simplify = F)
    Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
    Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]
    row.names(exp.mar) <- exp.mar$Protein
    exp.mar$Protein <- NULL
    exp.mar$Organism <- NULL
    design$sample_name <- str_remove(design$sample_name, fixed("_"))
    colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))
    
    exp.mar <- exp.mar[,colnames(exp.mar) %in% design$sample_name]
    design <- design[design$sample_name %in% colnames(exp.mar),]
    
    exp.mar <- exp.mar[,design$sample_name]
    
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
           df.temp <- log2(df.temp+1)
           df.temp <- seqKNNimp(df.temp)
           df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
       }else{
           df.temp[df.temp == 0 ] <- NA
           df.temp <- log2(df.temp+1)
           df.temp <- seqKNNimp(df.temp)
           df.temp <- data.frame(df.temp, check.names = F, check.rows = F)
        }

        

        sample_info <- data.frame(
            sample_id = colnames(df.temp),
            group = c(rep(i.Contrasts1, select.col.Contrasts1),
                      rep(i.Contrasts2, select.col.Contrasts2)))
            
        rownames(sample_info) <- row.names(sample_info$sample_id)
        sample_info$group <- as.factor(sample_info$group)
        
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
        meta.info <- c("group")
        niter <- 500 # Number of bootstrap samples
        K <- nrow(feature_info) / 2 # Set the value for K based on the data size
        K <- floor(K)
        group.name <- "group"
        formula.str <- "~ 0 + group" # Formula for group comparison
        
        # Run LimROTS analysis with trend and robust settings enabled
        se <- LimROTS(
            x = se,
            niter = niter, K = K, meta.info = meta.info,
            BPPARAM  = SnowParam(5, progressbar = T), group.name = group.name,
            formula.str = formula.str, trend = TRUE,
            robust = TRUE, permutating.group = FALSE
        )
        
        saveRDS(se, paste0("FragPipe_results/", "LimROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
     
        remove(se)
        gc()
        
        rots_results <- ROTS(data = df.temp , groups = as.numeric(sample_info$group), B = niter, K = K, 
                             seed = 1597, progress = TRUE, verbose = TRUE)
        
        saveRDS(rots_results, paste0("FragPipe_results/", "ROTS_" , i.exp , "_" , i.Contrasts, ".rds"))
        
        remove(rots_results)
        gc()
        
        design.model = model.matrix(~0+group , data = sample_info) 
        
        fit1 = lmFit(df.temp, design = design.model)
        cont <- makeContrasts(contrasts =  paste0( "group", i.Contrasts1 , "-" , "group" , i.Contrasts2), levels = design.model)
        fit2 = contrasts.fit(fit1,contrasts = cont)
        fit3 <- eBayes(fit2, trend = T, robust = T)
        limma.results = topTable(fit3, adjust="BH", sort.by = 'logFC', n=Inf)
        
        saveRDS(limma.results, paste0("FragPipe_results/", "Limma_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(limma.results)
        gc()
        
        
        
        
        # run SAM
        
        data<-list(x=as.matrix(df.temp), y=  sample_info$group
                   , geneid=feature_info$protein_id,
                   genenames=feature_info$protein_id, logged2=TRUE)
        
        
        samr.obj <- tryCatch({ samr( data, resp.type="Two class unpaired", nperms=100)
            }, error = function(e){
                return(NA)
            })
        
        if ( inherits(samr.obj, "list") ){
            
                             
        logFC<-log2(samr.obj$foldchange)
        
        pvalue <- data.frame(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
        
        adj.pvalue <- p.adjust(pvalue$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar., method = "BH")
        
        SAM.results<-cbind(logFC, pvalue, adj.pvalue)
        
        saveRDS(SAM.results, paste0("FragPipe_results/", "SAM_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(SAM.results)
        gc()
        
        }
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
        
        
        saveRDS(i.test.result, paste0("FragPipe_results/", "ttest_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(i.test.result)
        gc()
        
        # run ANOVA
        
        
        ## codes from:
        #https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#compare-p-values-from-deqms-to-ordinary-t-test-anova-and-limma
        design.model = model.matrix(~0+group , data = sample_info) 
        
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
        saveRDS(anova.results, paste0("FragPipe_results/", "ANOVA_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(anova.results)
        gc()
        
        }
    }
}




