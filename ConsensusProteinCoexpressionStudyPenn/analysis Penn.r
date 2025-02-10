library(stringr)
library(rrcovNA)
library(limma)
library(sva)
library(ggplot2)
library(limma)
library(SummarizedExperiment)
library(LimROTS)
library(BiocParallel)
library(ROTS)
library(dplyr)

df <- read.csv("1.MaxQuant-proteinGroups.txt syn21449513.1.txt", sep = "\t")

df <- df[df$Razor...unique.peptides >= 2,]

row.names(df) <- df$Protein.IDs

df <- df[,grepl("LFQ.intensity" , colnames(df) ,fixed = TRUE)]

colnames(df) <- str_remove(colnames(df) , "LFQ.intensity." )

df <- df[!grepl( "CON__" , row.names(df), fixed = TRUE),]
df <- df[!grepl( "REV__" , row.names(df), fixed = TRUE),]


meta <- read.csv("0.Traits Penn.csv")

table(meta$Dx.Final)

meta <- meta[meta$Dx.Final %in% c("Alzheimer's disease" , "Control"),]

meta <- meta[!duplicated(meta$Sample.ID),]

meta <- meta[order(meta$Dx.Final),]



colnames(df) <- str_split_fixed(colnames(df) , fixed("_"), 4 )[,3]
colnames(df) <- as.character( as.integer(colnames(df)) )
df <- df[,colnames(df) %in% meta$Sample.ID]
colnames(df) <- paste0("sample", colnames(df))
meta$Sample.ID <- paste0("sample", meta$Sample.ID)
meta <- meta[order(meta$Dx.Final),]
df <- df[,meta$Sample.ID]

df[df == 0] <- NA

df <- log2(df+1)

hist(as.matrix(df))

# Calculate the percentage of missing values for each row
row_missing_percent <- rowSums(is.na(df)) / ncol(df) * 100

# Remove rows with more than 50% missing values
df <- df[row_missing_percent <= 50, ]

row_missing_percent <- rowSums(is.na(df)) / ncol(df) * 100


df <- impSeq(df)

hist(as.matrix(df))

table(meta$Dx.Final)

meta <- meta[order(meta$Dx.Final),]
df <- df[,meta$Sample.ID]

plotMDS(df, labels = meta$Batch)

df.batched <- ComBat(df, batch = meta$Batch)

plotMDS(df.batched, labels = meta$Batch)

################################
#######
meta <- meta[order(meta$Dx.Final),]

df <- data.frame(df, check.rows = F, check.names = F)
df <- df[,meta$Sample.ID]

df.batched <- data.frame(df.batched, check.rows = F, check.names = F)
df.batched <- df.batched[,meta$Sample.ID]


###############  ############### ###############  ###############

meta.group <- meta[meta$Group %in% c("Control" , "AD"),]
df.batched.group <- df.batched[,colnames(df.batched) %in% meta.group$Sample.ID]



meta.group <- meta.group[,c(2,15,3,4,6,7,11)]

row.names(meta.group)<-  meta.group$Sample.ID
meta.group$Sample.ID <- NULL

meta.group$Group.num <- ifelse(meta.group$Group == "AD" , 1 , 2)
meta.group$Group.num <- as.factor(meta.group$Group.num)
meta.group$Batch <- as.factor(meta.group$Batch)
meta.group$Race <- as.factor(meta.group$Race)
meta.group$Age <- ifelse(meta.group$Age == "90+" , 91 , meta.group$Age)
meta.group$Age <- as.integer(meta.group$Age)
str(meta.group)



model.mat <- model.matrix(~0 + Group.num + Race + Sex + RunOrder + Batch + Age, data = meta.group)

df.batched.group <- df.batched.group[,row.names(model.mat)]
meta.group <- meta.group[row.names(model.mat),]

meta.group <- meta.group[order(meta.group$Group.num),]
df.batched.group <- data.frame(df.batched.group, check.rows = F, check.names = F)
df.batched.group <- df.batched.group[,row.names(meta.group)]



df <- df[,row.names(model.mat)]
meta.group <- meta.group[row.names(model.mat),]

meta.group <- meta.group[order(meta.group$Group.num),]
df <- data.frame(df, check.rows = F, check.names = F)
df <- df[,row.names(meta.group)]

###############  ############### ###############  ###############


feature_info <- data.frame(
    gene_id = row.names(df)
)
rownames(feature_info) <- feature_info$protein_id


# Create the SummarizedExperiment object
se <- SummarizedExperiment(
    assays = list(log.df = df),
    colData = meta.group,
    rowData = feature_info
)



# Set metadata and formula for LimROTS analysis
meta.info <- colnames(colData(se))
niter <- 1000 # Number of bootstrap samples
K <- nrow(feature_info) / 4 # Set the value for K based on the data size
K <- floor(K)
group.name <- c("Group.num")
formula.str <- "~ 0 + Group.num + Race + Sex + RunOrder + Age + Batch" # Formula for group comparison

set.seed(1234, kind = "default" , sample.kind = "default")


se <- LimROTS(
    x = se,
    niter = niter, K = K, meta.info = meta.info,
    BPPARAM  = SnowParam(10, progressbar = T), group.name = group.name,
    formula.str = formula.str, trend = TRUE,
    robust = TRUE, permutating.group = FALSE
)

saveRDS(se, paste0("PennResult//", "LimROTS_ADvscontrol.test" , ".rds"))
remove(se)
gc()



#### ROTS


set.seed(1234, kind = "default" , sample.kind = "default")


rots_results <- ROTS(data = df.batched.group , groups = as.numeric(meta.group$Group.num), 
                     B = niter, K = K, 
                     seed = 1234, progress = TRUE, verbose = TRUE)

saveRDS(rots_results, paste0("PennResult//", "ROTS_ADvsPD.batched.PDD" , ".rds"))

remove(rots_results)
gc()



model.mat <- model.matrix(~0 + Group.num + Race + Sex + RunOrder + Age + Batch, data = meta.group)
colnames(model.mat) <- make.names(colnames(model.mat))
fit1 <- lmFit(df.batched.group, model.mat)
cont <- makeContrasts(contrasts = c("Group.num1 - Group.num2"),
                      levels = model.mat)
fit2 = contrasts.fit(fit1,contrasts = cont)
fit3 <- eBayes(fit2, trend = T, robust = T)
limma.trend = topTable(fit3, adjust="BH" , n=Inf)
saveRDS(limma.trend, paste0("PennResult//", "limmatrend_ADvsPD.PDD" ,".rds"))
remove(limma.trend)
remove(fit1)
remove(fit2)
remove(model.mat)
remove(fit3)
gc()
