BiocManager::install("genomation")

setwd("/srv/data/my_shared_data_folder/epilepsy/data/bisulfiteseq")


# Main analysis package
library("methylKit")
# Annotation package
library("genomation")
library("GenomicRanges")

file.list=list("/srv/data/my_shared_data_folder/epilepsy/data/bisulfiteseq/SRR10493735_1.fa_bismark_bt2_pe.bismark.cov.gz", 
               "/srv/data/my_shared_data_folder/epilepsy/data/bisulfiteseq/SRR10493737_1.fa_bismark_bt2_pe.bismark.cov.gz")


## read the files to a methylRawList object: myobj
myobj=methRead(file.list,
               sample.id=list("SRR10493735",
                              "SRR10493737"),
               pipeline = "bismarkCoverage",
               assembly="hg38",
               treatment=c(1,0),
               context="CpG",
               mincov = 10
)

# check number of samples
myobj

# type of data stored here
head(myobj[[1]])

## Descriptive Statistics
getMethylationStats(myobj[[2]], plot=TRUE, both.strands=FALSE)

# Get a histogram of the read coverage per sample
getCoverageStats(myobj[[2]], plot=TRUE, both.strands=FALSE)

# Get percentile data by setting plot=FALSE
getCoverageStats(myobj[[1]], plot=FALSE, both.strands=FALSE)

## Filtration

myobj.filt <- filterByCoverage(myobj,
                               lo.count=10,
                               lo.perc=NULL,
                               hi.count=NULL,
                               hi.perc=99.9)

## Normalization

myobj.filt.norm <- normalizeCoverage(myobj.filt, method = "median")

## Merge Data

meth <- unite(myobj.filt.norm, destrand=FALSE)
#nrow(meth)/ 14399679

## Further Filtering

# get percent methylation matrix
pm=percMethylation(meth)

# calculate standard deviation of CpGs
sds=matrixStats::rowSds(pm)

# Visualize the distribution of the per-CpG standard deviation
# to determine a suitable cutoff
hist(sds, breaks = 100)

# keep only CpG with standard deviations larger than 2%
meth <- meth[sds > 2]

# This leaves us with this number of CpG sites
nrow(meth)
#10710104

## Remove known C -> T mutations

# give the locations of 2 example SNPs
mut <- GRanges(seqnames=c("chr21","chr21"),
               ranges=IRanges(start=c(9853296, 9853326),
                              end=c( 9853296,9853326)))

# select CpGs that do not overlap with mutations
meth <- meth[!as(meth,"GRanges") %over% mut, ]

## Data Structure/Outlier Detection/ check the correlation between samples
getCorrelation(meth,plot=TRUE)

# check the correlation between samples using dendrogram and PCA
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)

PCASamples(meth)

## Differential Methylation... This might take a few minutes.
myDiff <- calculateDiffMeth(meth,
                            overdispersion = "MN",
                            adjust="BH")
myDiff

# Simple volcano plot to get an overview of differential methylation
plot(myDiff$meth.diff, -log10(myDiff$qvalue))
abline(v=0)

# Overview of percentage hyper and hypo CpGs per chromosome.
diffMethPerChr(myDiff)

# get hyper methylated bases and order by qvalue
myDiff25p.hyper <- getMethylDiff(myDiff,
                                 difference=25,
                                 qvalue=0.01,
                                 type="hyper")
myDiff25p.hyper <- myDiff25p.hyper[order(myDiff25p.hyper$qvalue),]

# get hypo methylated bases and order by qvalue
myDiff25p.hypo <- getMethylDiff(myDiff,
                                difference=25,
                                qvalue=0.01,
                                type="hypo")
myDiff25p.hypo <- myDiff25p.hypo[order(myDiff25p.hypo$qvalue),]

# get all differentially methylated bases and order by qvalue
myDiff25p <- getMethylDiff(myDiff,
                           difference=25,
                           qvalue=0.01)
myDiff25p <- myDiff25p[order(myDiff25p$qvalue),]

## CpG Annotation

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)

files <- getSampleFiles()
print(files)

peak <- readPeakFile(files[[4]])
peak

covplot(peak, weightCol="V5")

covplot(peak, weightCol="V5", chrs=c("chr17", "chr18"), xlim=c(4.5e7, 5e7))

data("tagMatrixList")
tagMatrix <- tagMatrixList[[4]]

tagHeatmap(tagMatrix, xlim=c(-3000, 3000), color="red")

peakHeatmap(files[[4]], TxDb=txdb, upstream=3000, downstream=3000, color="red")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000),
            xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf2(files[[4]], TxDb=txdb, upstream=3000, downstream=3000,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")

plotAvgProf(tagMatrix, xlim=c(-3000, 3000), conf = 0.95, resample = 1000)

peakAnno <- annotatePeak(files[[4]], tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

plotAnnoPie(peakAnno)


plotAnnoBar(peakAnno)

upsetplot(peakAnno)

plotDistToTSS(peakAnno,
              title="Distribution of transcription factor-binding loci\nrelative to TSS")

library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)

data("tagMatrixList")
plotAvgProf(tagMatrixList, xlim=c(-3000, 3000))

plotAvgProf(tagMatrixList, xlim=c(-3000, 3000), conf=0.95,resample=500, facet="row")

tagHeatmap(tagMatrixList, xlim=c(-3000, 3000), color=NULL)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb,
                       tssRegion=c(-3000, 3000), verbose=FALSE)

plotAnnoBar(peakAnnoList)

plotDistToTSS(peakAnnoList)

genes = lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
names(genes) = sub("_", "\n", names(genes))
compKEGG <- compareCluster(geneCluster   = genes,
                           fun           = "enrichKEGG",
                           pvalueCutoff  = 0.05,
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 15, title = "KEGG Pathway Enrichment Analysis")

genes= lapply(peakAnnoList, function(i) as.data.frame(i)$geneId)
vennplot(genes)

p <- GRanges(seqnames=c("chr1", "chr3"),
             ranges=IRanges(start=c(1, 100), end=c(50, 130)))
shuffle(p, TxDb=txdb)

enrichPeakOverlap(queryPeak     = files[[5]],
                  targetPeak    = unlist(files[1:4]),
                  TxDb          = txdb,
                  pAdjustMethod = "BH",
                  nShuffle      = 50,
                  chainFile     = NULL,
                  verbose       = FALSE)

getGEOspecies()

getGEOgenomeVersion()

hg19 <- getGEOInfo(genome="hg19", simplify=TRUE)
head(hg19)

downloadGEObedFiles(genome="hg19", destDir="hg19")

gsm <- hg19$gsm[sample(nrow(hg19), 10)]
downloadGSMbedFiles(gsm, destDir="hg19")


