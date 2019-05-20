
args <- commandArgs(TRUE)

sessionInfo()
library('DESeq2')
library(ggplot2)
library(gplots)
library(tibble)

setwd(getwd())
#setwd('Y:/jupyter/UTR/new')
meta <- read.table(args[1], sep ='\t')
#meta <- read.table('metadat.txt', sep ='\t')
control = as.character(meta[1,1])
treat = as.character(meta[2,1])
rep1 = meta[1,2]
rep2 = meta[2,2]

data <- read.table(args[2], header = TRUE, sep ='\t')
#data <- read.table('UTRExtension_final.txt', header = TRUE, sep ='\t')

length <- length(colnames(data))

raw.data <- data[c(1, 8:length)]

newlength <- length(colnames(raw.data))
#head(raw.data)

rownames(raw.data) <- raw.data$UTRID

#dim(raw.data)

## half of the samples should have counts >=5
raw.data1 = raw.data[apply(raw.data[, -1], MARGIN = 1, function(x) {sum(x >= 5) >= length(x)/2}), ]

condition1 = factor(c(rep(control,rep1), rep(treat, rep2)))

countData1 <- as.matrix(raw.data1[,2:newlength], group = condition1)

#summary(countData1)

colData1 <- data.frame(row.names=colnames(raw.data1[,2:newlength]), group=condition1)

dds <- DESeqDataSetFromMatrix(countData=countData1, colData=colData1, design=~group)


dds$group <- factor(dds$group, levels = c(control, treat))
dds <- DESeq(dds)

res <- results(dds)

#res
res <- na.omit(res)
res1 <- res[res$padj<0.05,]
final <- res1[abs(res1$log2FoldChange) >=1,]

final <- as.data.frame(final)
#resultsNames(dds)


finalUTR <- rownames(final)

finalfile <- raw.data[finalUTR,]
rownames(finalfile) <- NULL

merged <- merge(data[1:7], finalfile, by.x = "UTRID", by.y = "UTRID")

final1 <- rownames_to_column(final, "UTRID")

finalMerged <- merge(merged,final1, by.x = "UTRID", by.y = "UTRID" )
write.table(finalMerged, args[3], sep = '\t', quote = FALSE, row.names = FALSE)
