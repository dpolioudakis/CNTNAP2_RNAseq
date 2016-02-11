# Damon Polioudakis
# 2016-02-05
# Normalize RNAseq data
#   (1) Counts + FPKM Data
#   (2) Filter Gene Expression
#   (3) View Data Pre-Normalzation and Pre-QC
#   (4) Normalize
#   (5) Remove Outliers
#   (6) Principal Component Analysis
#   (7) View + Filter Data Post-Normalization and Post-QC
# TODO
#   gcbias section is incomplete
################################################################################

rm(list = ls())
sessionInfo()

library(ggplot2)
library(reshape2)
library(xlsx)
library(DESeq2)
library(WGCNA)
library(cqn)

# Load data and assign variables

exDatDF <- read.csv("../data/HTSC/Exprs_HTSCexon.csv", row.names = 1)

seqQCdatDF <- read.csv("../metadata/PicardToolsQC.csv")

metDatDF <- read.xlsx("../metadata/Cntnap2 WT KO RNA-seq analysis.xlsx"
                      , sheetIndex = 1)
metDatDF <- metDatDF[ ,-1]

load("../source/Ensembl_Mm10_Exon_Union_Anno.rda")
mm10exUnionAnno <- mm10exUnionAnno

outGraphs <- "../analysis/graphs"
graphCodeTitle <- "Normalize_and_Outliers.R"
################################################################################

# Formatting and Filtering of Data

# Format expression data and metadata so that metadata rows and expression data
# columns are in same sample order
exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
metDatDF <- metDatDF[order(metDatDF$Sample.ID), ]

# Filter htseqcounts 

# Remove statistics from expression data frame
exDatDF <- exDatDF[-c((nrow(exDatDF)-4):nrow(exDatDF)), ]

# exp > 10 in 80% of samples
pres <- apply(exDatDF > 10, 1, sum)
idx <- which(pres > 0.8 * dim(exDatDF)[2])
fExDatDF <- exDatDF[idx, ]
################################################################################

# (3) View Data Pre-Normalization and Pre-QC

# Functions
# Function to output data frame of MDS PCs values and PCs
calcMDS <- function (exprDF) {
  # dist calculates the distances between the rows of a data matrix
  # Transpose data frame so samples are rows and genes are columns
  mds = cmdscale(dist(t(exprDF)), eig = T)
  pc1 = mds$eig[1]^2 / sum(mds$eig^2)
  pc2 = mds$eig[2]^2 / sum(mds$eig^2)
  mdsAndTreatmentLDF <- data.frame(mds$points
                                   , pc1 = pc1, pc2 = pc2)
  mdsAndTreatmentLDF
}

# MDS All Genes

# Calculate MDS
mdsDF <- calcMDS(fExDatDF)
# Add column with sample type info
mdsDF$type <- metDatDF$Genotype.factor

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 3) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1, size = 1.5) +
  scale_color_discrete(name = "Sample Type") +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq Counts", sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(paste0(outGraphs, "Normalize_and_Outliers_MDS_preN_E10in80.pdf"))

# Boxplot of counts for each sample
pdf(paste0(outGraphs, "Normalize_and_Outliers_BoxPlot_preN_E10in80.pdf"))
boxplot(log2(fExDatDF + 1), range = 0, col = as.factor(metDatDF$Genotype.factor)
        , main = paste(graphCodeTitle, "Boxplot HTSeq Counts", sep = "\n")
        , xlab = "", ylab = "Counts"
        , par(las = 2), par(mar=c(14,4,3,1)))
par(las = 1)
mtext("Sample", side = 1, line = 12)
dev.off()

# Density plot of counts for each sample
pdf(paste0(outGraphs, "Normalize_and_Outliers_Density_preN_E10in80.pdf"))
plot(density(log2(fExDatDF[ ,1] + 1))
     , col = as.factor(metDatDF$Genotype.factor)[1]
     , main = paste(graphCodeTitle, "Density HTSeq Counts", sep = "\n")
     , xlab = "Counts", ylim = c(0, 0.2))
for (i in 2:dim(fExDatDF)[2]) {
  lines(density(log2(fExDatDF[ ,i] + 1))
        , col = as.factor(metDatDF$Genotype.factor)[i])
}
dev.off()
################################################################################

# (4) Normalize

# Use CQN to normalize for GC content and gene length and quantile
RunCQN <- function (exprDatDF, gcLengthDF) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(gcLengthDF), rownames(exprDatDF))
  geneAnno <- gcLengthDF[match(keepV, rownames(gcLengthDF)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[ ,1] == 0] <- 1
  exprDatDF <- exprDatDF[match(keepV, rownames(exprDatDF)), ]
  
  # Run CQN with GC, length, and quantile normalization
  cqnDat <- cqn(exprDatDF, lengths = as.numeric(geneAnno[ ,1])
                , x = as.numeric(geneAnno[ ,2]), lengthMethod = c("smooth")
                , sqn = TRUE)
  # Get the log2(Normalized FPKM) values
  cqnDat <- cqnDat$y + cqnDat$offset
  cqnDat
}
cqnFexDat <- RunCQN(fExDatDF, mm10exUnionAnno)

# Check correlation of expression level to GC content and gene length pre and
# post CQN normalization

CheckCQNnorm <- function (preCQN, postCQN, gcLengthDF) {
  # Remove genes not in GC and gene length annotation file
  keepV <- intersect(rownames(gcLengthDF), rownames(preCQN))
  geneAnno <- gcLengthDF[match(keepV, rownames(gcLengthDF)), ]
  # Set genes with length 0 to length 1 - why? - from Vivek's code
  geneAnno[geneAnno[,1] == 0] <- 1
  keepgenes <- intersect(rownames(preCQN),rownames(postCQN))
  preCQN <- preCQN[match(keepgenes,rownames(preCQN)),]
  postCQN <- postCQN[match(keepgenes,rownames(postCQN)),]
  geneAnno1 <- geneAnno[match(keepgenes,rownames(geneAnno)),]
  
  qcCorrCQNm <- matrix(NA,nrow=ncol(preCQN),ncol=4)
  colnames(qcCorrCQNm) <- c("preNormGCcorr", "preNormLengthCorr"
                            ,"postNormGCcorr", "postNormLengthCorr")
  for (i in 1:nrow(qcCorrCQNm)) {
    qcCorrCQNm[i,1] <- cor(preCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,2] <- cor(preCQN[,i], geneAnno1[,1], method="spearman")
    qcCorrCQNm[i,3] <- cor(postCQN[,i], geneAnno1[,2], method="spearman")
    qcCorrCQNm[i,4] <- cor(postCQN[,i], geneAnno1[,1], method="spearman")
  }
  qcCorrCQNm
}

qcCorrCQNm <- CheckCQNnorm(exDatDF, cqnFexDat, mm10exUnionAnno)
apply(qcCorrCQNm, 2, quantile)
qcCorrCQNm <- data.frame(qcCorrCQNm)
qcCorrCQNm <- melt(qcCorrCQNm)
colnames(qcCorrCQNm) <- c("CorrType", "Corr")

# Histogram of spearman's rho pre and post normalization for GC and gene length
ggplot(qcCorrCQNm, aes(x = Corr)) +
  facet_wrap(~CorrType, nrow = 2) +
  geom_histogram(binwidth = 0.01) +
  ylab("Counts") +
  xlab("Spearman's rho across samples") +
  labs(title = paste(graphCodeTitle
                , "Histogram: Spearman's rho across samples - pre and post CQN"
                , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(file = paste0(outGraphs, "/Normalize_and_Outliers_PrePostCQN_E10in80.pdf")
       , height = 6)
################################################################################

# (5) Remove Outliers

# We will go with the less rigorous 'All Samples' method here, but the 'By Dx' method is also included for your reference

## For All Samples

pdf("./SampleNetwork_allSamps.pdf")


normadj <- (0.5+0.5*bicor(fExDatDF))^2 ## Calculate connectivity
normadj <- (0.5+0.5*bicor(cqnFexDat))^2 ## Calculate connectivity

netsummary <- fundamentalNetworkConcepts(normadj)
par(mfrow=c(2,1))
C <- netsummary$ClusterCoef
z.C <- (C-mean(C))/sqrt(var(C))
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))

datLabel <- paste(metDatDF$Genotype.factor)        ## You can add more datMeta information here, ie. sex or race, if wanted
plot(1:length(z.ku),z.ku)
text(1:length(z.ku),z.ku,label=datLabel,pos=4,cex=0.6)
abline(h= -2, col="red")
flag <- z.ku < -2

plot(z.ku,z.C)
text(z.ku,z.C,label=datLabel,pos=4,cex=0.6)

dev.off()

datMeta = datMeta[!flag,]
datExpr.htg <- datExpr.htg[,!flag]
datExpr.cuff <- datExpr.cuff[,!flag]








## raw statistics

ggDF <- data.frame(t(rbind(Genotype.factor = metDatDF$Genotype.factor, fExDatDF)))
ggDF$Sample <- row.names(ggDF)
ggDF[1:3,1:4]
ggDF <- melt(ggDF, id.vars = c("Genotype.factor", "Sample"))
head(ggDF)
colnames(ggDF) <- c("Genotype", "Sample", "Gene", "Count")

# compute lower and upper whiskers
ylim1 <- boxplot.stats(ggDF$Count)$stats[c(1, 5)]

# scale y limits based on ylim1

ggplot(ggDF, aes(x = Sample, y = Count)) +
  geom_boxplot(aes(col = as.factor("Genotype"))) +
  coord_cartesian(ylim = ylim1*1.1)