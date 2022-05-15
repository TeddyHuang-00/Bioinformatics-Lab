if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("DESeq2")
BiocManager::install("pasilla")

browseVignettes(package = "DESeq2")

### Input
# Count matrix input
library("pasilla")
pasCts <- system.file("extdata",
  "pasilla_gene_counts.tsv",
  package = "pasilla", mustWork = TRUE
)
pasAnno <- system.file("extdata",
  "pasilla_sample_annotation.csv",
  package = "pasilla", mustWork = TRUE
)
cts <- as.matrix(read.csv(pasCts, sep = "\t", row.names = "gene_id"))
coldata <- read.csv(pasAnno, row.names = 1)
coldata <- coldata[, c("condition", "type")]
coldata$condition <- factor(coldata$condition)
coldata$type <- factor(coldata$type)
# same name
rownames(coldata) <- sub("fb", "", rownames(coldata))
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
# same order
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
# construct a DESeqDataSet
library("DESeq2")
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData = coldata,
  design = ~condition
)
dds
ddsMF <- dds

# htseq-count input
directory <- system.file("extdata",
  package = "pasilla",
  mustWork = TRUE
)
sampleFiles <- grep("treated", list.files(directory), value = TRUE)
sampleCondition <- sub("(.*treated).*", "\\1", sampleFiles)
sampleTable <- data.frame(
  sampleName = sampleFiles,
  fileName = sampleFiles,
  condition = sampleCondition
)
sampleTable$condition <- factor(sampleTable$condition)
ddsHTSeq <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~condition
)
ddsHTSeq

### pre-filtering
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]
dds$condition <- factor(dds$condition, levels = c("untreated", "treated"))

### Differential expression analysis
dds <- DESeq(dds)
res <- results(dds)
res
res <- results(dds, contrast = c("condition", "treated", "untreated"))
res
mcols(res)$description
### p-values
resOrdered <- res[order(res$pvalue), ]
summary(res)
sum(res$padj < 0.1, na.rm = TRUE)
res05 <- results(dds, alpha = 0.05) # an adjusted p value below a given FDR cutoff, alpha
summary(res05)

### plot
library(EnhancedVolcano)
EnhancedVolcano(res,
  lab = rownames(res),
  x = "log2FoldChange",
  y = "pvalue",
  title = "treated versus untreated",
  pCutoff = 10e-32,
  FCcutoff = 0.5,
  pointSize = 3.0,
  labSize = 6.0
)

plotCounts(dds, gene = "FBgn0039155", intgroup = "condition")

### export results
write.csv(as.data.frame(resOrdered),
  file = "condition_treated_results.csv"
)


### multi-factor
colData(dds)
levels(ddsMF$type)
# change the levels of type so it only contains letters
levels(ddsMF$type) <- sub("-.*", "", levels(ddsMF$type))
levels(ddsMF$type)
# As condition is the variable of interest, we put it at the end of the formula.
# Thus the results function will by default pull the condition results unless contrast or name arguments are specified.
design(ddsMF) <- formula(~ type + condition)
ddsMF <- DESeq(ddsMF)
resMF <- results(ddsMF)
head(resMF)
resMFType <- results(ddsMF,
  contrast = c("type", "single", "paired")
)
head(resMFType)