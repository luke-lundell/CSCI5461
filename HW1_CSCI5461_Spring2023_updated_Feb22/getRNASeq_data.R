#Deseq2 pipeline/TCGA data from:
#https://www.bioconductor.org/packages/release/data/experiment/vignettes/GSE62944/inst/doc/GSE62944.html#case-study
#for GBM:
library(ExperimentHub)
library(DESeq2)

eh = ExperimentHub()
query(eh , "GSE62944")
tcga_data <- eh[["EH1"]]

gbm_data <- tcga_data[, which(phenoData(tcga_data)$CancerType=="GBM")]


short_survival_idx <- which(as.numeric(phenoData(gbm_data)$death_days_to) <= 365)
short_survival_data <- exprs(gbm_data)[, short_survival_idx]

long_survival_idx <- which(as.numeric(phenoData(gbm_data)$death_days_to) > 365)
long_survival_data <- exprs(gbm_data)[, long_survival_idx]


# make a countTable.
countData <- cbind(short_survival_data, long_survival_data)
survival_times= c(as.numeric(phenoData(gbm_data)$death_days_to)[short_survival_idx],as.numeric(phenoData(gbm_data)$death_days_to)[long_survival_idx])

# for DE analysis with DESeq2 we need a sampleTable
samples= c(colnames(short_survival_data), colnames(long_survival_data))
group =c(rep("short",length(short_survival_idx)), rep("long", length(long_survival_idx)))
coldata <- cbind(samples, group)
colnames(coldata) <- c("sampleName", "Group")
coldata[,"Group"] <- factor(coldata[,"Group"], c("short","long"))
coldata <- DataFrame(coldata)

# Now we can run DE analysis

ddsMat <- DESeqDataSetFromMatrix(countData = countData,
                                 colData = coldata,
                                 design = ~ Group)

dds <- ddsMat
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds <- DESeq(dds)


res <- results(dds) 
summary(res)

coldata$Group <- factor(coldata$Group , levels=c("1", "2"), labels=c("short", "long"))

res <- res[ , !(names(res) %in% c("padj"))]

write.csv(countData, file = "HW1-GSE62944-count.csv")
write.csv(coldata, file = "HW1-GSE62944-clinical.csv")
write.csv(res, file = "HW1-DESeq2.csv")
