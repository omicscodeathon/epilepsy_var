
################################################################
#  Part 1 : Differential Genes Expression Analysis
################################################################

#Install required packages

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install('EnhancedVolcano')
BiocManager::install("apeglm")
BiocManager::install("ReportingTools")
BiocManager::install("goseq")
BiocManager::install("vsn")
install.packages("tidyverse")
install.packages("genefilter")
install.packages("pheatmap")
install.packages("RColorBrewer")
install.packages("PoiClaClu")
install.packages("ggplot2")
install.packages("ggbeeswarm")




### Load required R packages

library(DESeq2)
library(genefilter)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(RColorBrewer)
library(PoiClaClu)
library(ggplot2)
library(ggrepel)
library(EnhancedVolcano)
library(dplyr)
library(ggbeeswarm)
library(apeglm)
library(ReportingTools)
library(vsn)

##########################################################
##     STEP 1: preparing reads or fragment counts       ##


### Load features counts matrix dataset 
counts_data = read.table("featureCounts.txt", header = TRUE, sep = ",", stringsAsFactors = FALSE)

counts = counts_data

# loadind sample info
colData = read.delim("rnaseq_metadata.csv", sep = ",", stringsAsFactors = TRUE)

colnames(counts)[7:42] <- substr(colnames(counts)[7:42],start=51,stop=57)
counts = counts[,c(1,7:42)]
rownames(counts) =  counts$Geneid
counts = counts[,-c(1)]

#colData = colData[, c(4:6)]
rownames(colData) =  colData$Accession
colData = colData[,-c(1)]


colnames(counts) = rownames(colData)

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts) %in% rownames(colData))

# are they in the same order?
all(colnames(counts) == rownames(colData))

############################################################
##   STEP:2 Creating DESeq object from counts Matrix      ##


## -----------------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                                 colData = colData,
                                 design = ~ Tissue + Condition)

## --Pre-filtering the dataset--------------------------------------------------
nrow(dds)
keep <- rowSums(counts(dds)) > 10
dds <- dds[keep,]
nrow(dds)

###############################################################################
# STEP 3 using DESeq2 for differential expression analysis


## DESeq2 differential expression with DESeq function ##
dds <- DESeq(dds)

# extract results from dds DESeq2 results function

res <- results(dds, contrast=c("Tissue","Hippocampus","Neocortex"), alpha = 0.05)

## # summary of the results
summary(res)
head(res)


# DESeq2 performs for each gene a hypothesis test to see whether evidence 
# is sufficient to decide against the null hypothesis that there is zero effect 
# of the treatment on the gene and that the observed difference between 
# treatment and control was merely caused by experimental variability
resultsNames(dds)

resLFC <- lfcShrink(dds, coef="Tissue_Neocortex_vs_Hippocampus", type="apeglm", lfcThreshold=1)


## -----------------------------------------------------------------------------
mcols(res, use.names = TRUE)

####### MA plots? M (log ratio) and A (mean average)  ####### 
# overview of the distribution of the estimated coefficients, or comparisons of interest, across all genes

plotMA(res, ylim=c(-2,2))
plotMA(resLFC, ylim=c(-2,2))


# counts plot of individual genes
library("ggbeeswarm")

geneCounts <- plotCounts(dds, gene=which.min(res$padj), intgroup=c("Tissue","Condition"), returnData = TRUE)

plotCounts(dds, gene=which.min(res$padj), intgroup="Tissue")

countplt <- ggplot(geneCounts, aes(x=Tissue, y=count, color = Condition, group = Condition)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

countplt


# extract VST transformed value for all genes for each sample
library("vsn")
vsd <- vst(dds, blind=FALSE)

# generate heatmap

select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]


df <- as.data.frame(colData(dds)[,c("Tissue","Condition")])

# heatmap of top 20 significant genes

pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=TRUE,
         cluster_cols=TRUE, annotation_col=df)


## Heatmap of the sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))

library("RColorBrewer")
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$Tissue, vsd$Condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)

# Significance ------------------------------------------------------------
## subset results, consider a fraction of 10% false positives acceptable FDR alpha=0.1

resSig <- subset(res, res$padj < 0.1 & res$baseMean > 50)
summary(resSig)


# extract up regulated genes
resSigUp <- subset(resSig, log2FoldChange > 0)
resSigUp[(order(resSigUp$pvalue)) ,]


# extract down regulated
resSigDown <- subset(resSig, log2FoldChange < 0)
resSigDown

resOrdered <- resSig[(order(resSig$pvalue, decreasing = TRUE)) ,]
head(resOrdered)

## ---Annotating and exporting results------------------------------------

columns(org.Hs.eg.db)

## -----------------------------------------------------------------------------
resOrdered$symbol <- mapIds(org.Hs.eg.db,
                        keys=rownames(resOrdered),
                        column="SYMBOL",
                        keytype="ENSEMBL",
                        multiVals="first")
resOrdered$entrez <- mapIds(org.Hs.eg.db,
                        keys=rownames(resOrdered),
                        column="ENTREZID",
                        keytype="ENSEMBL",
                        multiVals="first")
head(resOrdered)

## -----volcano plot-------------------------------------------------------

EnhancedVolcano(resOrdered,
                lab = resOrdered$symbol,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Hippocampus versus Neocortex',
                labSize = 6.0)

## ----Exporting results--------------------------------------------------------
resOrderedDF <- as.data.frame(resOrdered) 

write.csv(resOrderedDF, file = "DE_results.csv")


htmlRep <- HTMLReport(shortName="report", title="My report",
                       reportDirectory="./report")
publish(resOrderedDF, htmlRep)
url <- finish(htmlRep)
browseURL(url)

#############################################################
# Part 2 : Functional enrichment analysis of the DE genes
#############################################################
BiocManager::install("clusterProfiler")
library(goseq)
library(clusterProfiler)


gene_to_test = resOrdered[resOrdered$log2FoldChange > 0.5,]
gene_to_test = gene_to_test$symbol

go_resuts_PB = enrichGO(gene = gene_to_test, OrgDb = "org.Hs.eg.db", keyType = "SYMBOL",
                     ont = "BP")

go_plot = plot(barplot(go_resuts_PB, showCategory = 20))

go_plot

png("go_enrich.png", res = 250, width = 1200, height = 1800)
print(go_plot)
dev.off()

go_df <- as.data.frame(go_resuts_PB) 

write.csv(go_df, file = "GO_results.csv")
