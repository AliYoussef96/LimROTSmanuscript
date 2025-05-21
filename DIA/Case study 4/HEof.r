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

data(UPS1.Case4)

colnames(colData(UPS1.Case4))

set.seed(1597, sample.kind = "default" , kind = "default")

        
        # Set metadata and formula for LimROTS analysis
        meta.info <- c("tool" , "Conc.", "fake.batch")
        niter <- 1000 # Number of bootstrap samples
        K <- nrow(assay(UPS1.Case4)) / 4 # Set the value for K based on the data size
        K <- floor(K)
        group.name <- "Conc."
        formula.str <- "~0+Conc.+tool+fake.batch" # Formula for group comparison
        
        # Run LimROTS analysis with trend and robust settings enabled
        UPS1.Case4 <- LimROTS(
            x = UPS1.Case4,
            niter = niter, K = K, meta.info = meta.info,
            BPPARAM  = SnowParam(10, progressbar = T), group.name = group.name,
            formula.str = formula.str, trend = TRUE,
            robust = TRUE, permutating.group = FALSE
        )
        
        saveRDS(UPS1.Case4, paste0("results/", "LimROTS" ,  ".rds"))
        
        remove(UPS1.Case4)
        gc()
        
        data(UPS1.Case4)
        
        rots_results <- ROTS(data = assay(UPS1.Case4) , groups = as.numeric(UPS1.Case4$Conc.), B = niter, K = K, 
                             seed = 1597, progress = TRUE, verbose = TRUE)
        
        saveRDS(rots_results, paste0("results/", "ROTS" ,  ".rds"))
        
        remove(rots_results)
        gc()
        
        
        rots_results <- lmROTS(formula = "Conc.+tool+fake.batch",data = assay(UPS1.Case4) , 
                               metadata  = data.frame(colData(UPS1.Case4)), B = niter, K = K, 
                             seed = 1597 , BPPARAM  = SnowParam(10, progressbar = T))
        
        saveRDS(rots_results, paste0("results/", "ROTS.lm" ,  ".rds"))
        
        remove(rots_results)
        gc()
        
        sample_info <- colData(UPS1.Case4)
        
        design.model = model.matrix(~0+Conc.+tool+fake.batch , data = sample_info) 
        colnames(design.model) <- make.names(colnames(design.model))
        fit1 = lmFit(assay(UPS1.Case4), design = design.model)
        cont <- makeContrasts(contrasts =  "Conc.2-Conc.1", levels = design.model)
        fit2 = contrasts.fit(fit1,contrasts = cont)
        fit3 <- eBayes(fit2, trend = T, robust = T)
        limma.results = topTable(fit3, adjust="BH", coef = "Conc.2-Conc.1"  ,sort.by = 'logFC', n=Inf)
        
        saveRDS(limma.results, paste0("results/", "Limma" ,  ".rds"))
        remove(limma.results)
        gc()
        
        feature_info <- rowData(UPS1.Case4)
        
        # run SAM
        
        data<-list(x=as.matrix(assay(UPS1.Case4)), y=sample_info$Conc.
                   , geneid=feature_info$GeneID,
                   genenames=feature_info$GeneID, logged2=TRUE)
        
        samr.obj <- samr( data, resp.type="Two class unpaired", nperms=100)
        logFC<-log2(samr.obj$foldchange)
        
        pvalue <- data.frame(samr.pvalues.from.perms(samr.obj$tt, samr.obj$ttstar))
        
        adj.pvalue <- p.adjust(pvalue$samr.pvalues.from.perms.samr.obj.tt..samr.obj.ttstar., method = "BH")
        
        SAM.results<-cbind(logFC, pvalue, adj.pvalue)
        
        saveRDS(SAM.results, paste0("results/", "SAM" ,  ".rds"))
        remove(SAM.results)
        gc()
        
        df.temp <- assay(UPS1.Case4)
        # run t-test
        i.test.result <- data.frame()
        for(i.test in seq_len(nrow(df.temp))){
            
            i.test.df <- df.temp[i.test,]
            
            
            pvalue <- tryCatch({
                t.test(formula = as.numeric(i.test.df)~sample_info$Conc.)$p.value},
                error = function(e){
                    1
                })
            
            #fc.calc <- mean( as.numeric( i.test.df[,1:select.col.Contrasts1] ) ) - mean( as.numeric(i.test.df[,select.col.Contrasts1+1:select.col.Contrasts2]) )
            fc.calc <- NA
            
            i.test.result <- rbind(i.test.result, data.frame(row.names = row.names(i.test.df), 
                                                             pvalue = pvalue , 
                                                             fc.calc = fc.calc))
        }
        
        i.test.result$adj.pvalue <- p.adjust(i.test.result$pvalue, method = "BH")
        
        
        saveRDS(i.test.result, paste0("results/", "ttest" ,  ".rds"))
        remove(i.test.result)
        gc()
        
        # run ANOVA
        
        
        ## codes from:
        #https://www.bioconductor.org/packages/release/bioc/vignettes/DEqMS/inst/doc/DEqMS-package-vignette.html#compare-p-values-from-deqms-to-ordinary-t-test-anova-and-limma
        design.model = model.matrix(~0+Conc.+tool+fake.batch , data = sample_info) 
        
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
        saveRDS(anova.results, paste0("results/", "ANOVA" ,  ".rds"))
        remove(anova.results)
        gc()
        
        





