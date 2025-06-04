
# Load required packages
library(stringr)
library(SummarizedExperiment)
library(LimROTS)
library(ROTS)
library(BiocParallel)
library(rrcovNA)
library(imputeLCMD)
library(limma)
library(samr)

# List all experiments from the Spectronaut directory
all.exp <- list.files("Spectronaut/")
all.exp <- str_split_fixed(all.exp, "_DIA_" , 2)[,1]
all.exp <- unique(all.exp)

# Set seed for reproducibility
set.seed(1597, sample.kind = "default", kind = "default")

# Loop through each experiment folder
for(i.exp in all.exp){
  message("Processing experiment: ", i.exp)
  
  # Load intensity and design tables
  exp.mar <- read.csv(file.path("Spectronaut", paste0(i.exp, "_LFQ_Spectronaut_dlfq_pro_intensity.tsv")), sep = "\t")
  design <- read.csv(file.path("Spectronaut", paste0(i.exp, "_LFQ_Spectronaut_design.tsv")), sep = "\t")

  # Create condition contrasts (e.g., AB, AC)
  Contrasts <- combn(design$condition, 2, simplify = FALSE)
  Contrasts <- unique(sapply(Contrasts, function(x) paste(sort(x), collapse = "")))
  Contrasts <- Contrasts[!sapply(Contrasts, function(x) length(unique(strsplit(x, "")[[1]])) == 1)]

  # Preprocess expression data
  row.names(exp.mar) <- exp.mar$Protein
  exp.mar$Protein <- NULL
  exp.mar$Organism <- NULL
  design$sample_name <- str_remove(design$sample_name, fixed("_"))
  colnames(exp.mar) <- str_remove(colnames(exp.mar), fixed("_"))

  exp.mar <- exp.mar[, colnames(exp.mar) %in% design$sample_name]
  design <- design[design$sample_name %in% colnames(exp.mar),]
  exp.mar <- exp.mar[, design$sample_name]

  # Loop through each contrast
  for(i.Contrasts in Contrasts){
    message("Running contrast: ", i.Contrasts)

    # Define groups
    groups <- str_split(i.Contrasts, "")[[1]]
    i.Contrasts1 <- groups[1]
    i.Contrasts2 <- groups[2]

    select.col.Contrasts1 <- sum(design$condition == i.Contrasts1)
    select.col.Contrasts2 <- sum(design$condition == i.Contrasts2)

    # Require at least two samples per group
    if(select.col.Contrasts1 > 1 & select.col.Contrasts2 > 1){
      design.temp <- design[design$condition %in% groups,]
      df.temp <- exp.mar[, design.temp$sample_name]

      # Log2 transform and impute missing
      df.temp[df.temp == 0] <- NA
      df.temp <- log2(df.temp + 1)
      df.temp <- impSeq(df.temp)
      df.temp <- data.frame(df.temp, check.names = FALSE)

      # Create sample and feature info
      sample_info <- data.frame(
        sample_id = colnames(df.temp),
        group = factor(c(rep(i.Contrasts1, select.col.Contrasts1),
                         rep(i.Contrasts2, select.col.Contrasts2)))
      )
      rownames(sample_info) <- sample_info$sample_id

      feature_info <- data.frame(
        protein_id = rownames(df.temp)
      )
      rownames(feature_info) <- feature_info$protein_id

      # Create SummarizedExperiment object
      se <- SummarizedExperiment(
        assays = list(protein.exp = df.temp),
        colData = sample_info,
        rowData = feature_info 
      )

      # Run LimROTS
      se <- LimROTS(
        x = se,
        niter = 1000,
        K = floor(nrow(df.temp) / 2),
        meta.info = "group",
        BPPARAM = SnowParam(5, progressbar = TRUE),
        group.name = "group",
        formula.str = "~ 0 + group",
        trend = TRUE,
        robust = TRUE
      )
      saveRDS(se, file.path("Spectronaut_results", paste0("LimROTS_", i.exp, "_", i.Contrasts, ".rds")))
      rm(se); gc()

      # ROTS
      rots_results <- ROTS(
        data = df.temp,
        groups = as.numeric(sample_info$group),
        B = 1000,
        K = floor(nrow(df.temp) / 2),
        seed = 1597,
        progress = TRUE,
        verbose = TRUE
      )
      saveRDS(rots_results, file.path("Spectronaut_results", paste0("ROTS_", i.exp, "_", i.Contrasts, ".rds")))
      rm(rots_results); gc()

      # Limma
      design.model <- model.matrix(~0 + group, data = sample_info)
      fit <- lmFit(df.temp, design = design.model)
      contrast <- makeContrasts(contrasts = paste0("group", i.Contrasts1, "-group", i.Contrasts2), levels = design.model)
      fit <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
      limma.res <- topTable(fit, adjust = "BH", sort.by = "logFC", n = Inf)
      saveRDS(limma.res, file.path("Spectronaut_results", paste0("Limma_", i.exp, "_", i.Contrasts, ".rds")))
      rm(limma.res); gc()

      # SAM
      sam.data <- list(
        x = as.matrix(df.temp),
        y = as.numeric(sample_info$group),
        geneid = feature_info$protein_id,
        genenames = feature_info$protein_id,
        logged2 = TRUE
      )
      sam.obj <- samr(sam.data, resp.type = "Two class unpaired", nperms = 100)
      logFC <- log2(sam.obj$foldchange)
      pval <- samr.pvalues.from.perms(sam.obj$tt, sam.obj$ttstar)
      adj.p <- p.adjust(pval, method = "BH")
      SAM.res <- cbind(logFC, pvalue = pval, adj.pvalue = adj.p)
      saveRDS(SAM.res, file.path("Spectronaut_results", paste0("SAM_", i.exp, "_", i.Contrasts, ".rds")))
      rm(SAM.res); gc()

      # t-test
      ttest.res <- t(sapply(1:nrow(df.temp), function(i){
        test <- t.test(as.numeric(df.temp[i,]) ~ sample_info$group)
        fc <- mean(df.temp[i, sample_info$group == i.Contrasts1]) - mean(df.temp[i, sample_info$group == i.Contrasts2])
        c(pval = test$p.value, fc = fc)
      }))
      ttest.res <- as.data.frame(ttest.res)
      ttest.res$adj.pvalue <- p.adjust(ttest.res$pval, method = "BH")
      rownames(ttest.res) <- rownames(df.temp)
      saveRDS(ttest.res, file.path("Spectronaut_results", paste0("ttest_", i.exp, "_", i.Contrasts, ".rds")))
      rm(ttest.res); gc()

      # ANOVA
        anova.results <- data.frame()
        for(i.test in seq_len(nrow(df.temp))){
          i.test.df <- df.temp[i.test,]
          fit = aov(as.numeric(i.test.df)~group, data = sample_info)
          fit <- summary(fit)[[1]]
          fc.calc <- mean( as.numeric( i.test.df[,1:select.col.Contrasts1] ) ) - mean( as.numeric(i.test.df[,select.col.Contrasts1+1:select.col.Contrasts2]) )
          anova.results = rbind(anova.results, data.frame(row.names = row.names(i.test.df),
                                                          logFC=log2FC,
                                                          P.Value=fit$`Pr(>F)`[1]) ) }
        anova.results$adj.P.Val <- p.adjust(anova.results$pvalue, method  = "BH")
        saveRDS(anova.results, paste0("spectronaut_resutls/", "ANOVA_" , i.exp , "_" , i.Contrasts, ".rds"))
        remove(anova.results)
        gc()
    }
  }
}
