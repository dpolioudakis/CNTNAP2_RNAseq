# Damon Polioudakis
# 2016-02-05
# Normalize RNAseq data
#   (1) Counts + FPKM Data
#   (2) Filter Gene Expression
#   (3) View Data Pre-Normalzation and Pre-QC
#   (4) Normalize
#   (5) Remove Outliers
#   (6) Principal Component Analysis
#   (7) View Data Post-Normalization and Post-QC
#   Output filtered, normalized, and outliers removed expression values
#   cqnFexDat to ../data/Expr_E10in80_Nm_OR.rda
# TODO
#   gcbias section is incomplete
#   Add Covariates PCA
#   Understand outlier removal
#   Save processed data
################################################################################

rm(list = ls())
sessionInfo()

library(ggplot2)
library(reshape2)
library(xlsx)
library(DESeq2)
library(WGCNA)
library(cqn)
library(corrplot)
library(limma)

# Load data and assign variables

exDatDF <- read.csv("../data/HTSC/Exprs_HTSCexon.csv", row.names = 1)

# Picard QC data
seqQCdatDF <- read.csv("../metadata/PicardToolsQC.csv")

metDatDF <- read.xlsx("../metadata/Cntnap2 WT KO RNA-seq analysis.xlsx"
                      , sheetIndex = 1)
metDatDF <- metDatDF[ ,-1]

# Gene Length and GC bias
load("../analysis/tables/Gene_Length_GC_Bias.rda")
lenBiasDF <- data.frame(GeneLengthScore = lenBias)
gcBiasDF <- data.frame(GCcontent = gcBias)

load("../source/Ensembl_Mm10_Exon_Union_Anno.rda")
mm10exUnionAnno <- mm10exUnionAnno

outGraphs <- "../analysis/graphs"
graphCodeTitle <- "Normalize_and_Outliers.R"
################################################################################

### Functions

## Functions for MDS analysis

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

## Functions for PC Analysis of Technical Covariates

# Run PCA funtion
Run_PCA <- function (data, pClabel = "Unregressed") {
  pCdat <- prcomp(data, center=F);
  # pCdat <- prcomp(exDatDF, center=F);
  topPCs <- pCdat$rotation[,1:5];
  # Calculate variance explained by each PC
  varExp <- (pCdat$sdev)^2 / sum(pCdat$sdev^2)
  topVar <- varExp[1:5]
  colnames(topPCs) <- paste(pClabel, "\n", colnames(topPCs)
                            , " (", signif(100 * topVar[1:5], 2), "%)", sep = "")
  return(topPCs)
}

# Correlation panel function
panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  if (class(x) == "numeric" & class(y) == "numeric") {
    r <- abs(cor(x, y, use = "pairwise.complete.obs", method = "spearman"))
  } else {
    lmout <- lm(y~x)
    r <- sqrt(summary(lmout)$adj.r.squared)
  }
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 2/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor * (r/3))
}

# Plot correlation panel function
Plot_Corr_Matrix <- function (data, color, title, outPath) {
  pdf(outPath, height = 20, width = 24)
  pairs(data, pch = 19, col = color
        , upper.panel = panel.cor
        , main = paste("Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values",
                       "\n", title, sep = "")
        # , cex.labels = (cex = 1.25)
  )
  dev.off()
}
################################################################################

### Formatting and Filtering of Data

## Formatting

# Add length bias to metadata
row.names(lenBiasDF) <- substring(row.names(lenBiasDF), 2)
metDatDF <- merge(x = metDatDF, y = lenBiasDF
                  , by.x = "Sample.ID", by.y = "row.names")
# Add GC bias to metadata
row.names(gcBiasDF) <- substring(row.names(gcBiasDF), 2)
metDatDF <- merge(x = metDatDF, y = gcBiasDF
                  , by.x = "Sample.ID", by.y = "row.names")

# Format expression data and metadata so that metadata rows and expression data
# columns are in same sample order
exDatDF <- exDatDF[ ,order(colnames(exDatDF))]
metDatDF <- metDatDF[order(metDatDF$Sample.ID), ]

# Remove statistics from expression data frame
exDatDF <- exDatDF[-c((nrow(exDatDF)-4):nrow(exDatDF)), ]

## Filter htseqcounts 

# exp > 10 in 80% of samples
pres <- apply(exDatDF > 10, 1, sum)
idx <- which(pres > 0.8 * dim(exDatDF)[2])
fExDatDF <- exDatDF[idx, ]
################################################################################

### (3) View Data Pre-Normalization and Pre-QC

## MDS All Genes

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
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq log2(Counts + 1)"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(paste0(outGraphs, "/Normalize_and_Outliers_MDS_preN_E10in80.pdf"))

######

## Boxplot of counts for each sample

pdf(paste0(outGraphs, "/Normalize_and_Outliers_BoxPlot_preN_E10in80.pdf"))
boxplot(log2(fExDatDF + 1), range = 0, col = as.factor(metDatDF$Genotype.factor)
        , main = paste(graphCodeTitle, "Boxplot HTSeq log2(Counts + 1)"
                       , sep = "\n")
        , xlab = "", ylab = "Counts"
        , par(las = 2), par(mar=c(14,4,3,1)))
par(las = 1)
mtext("Sample", side = 1, line = 12)
dev.off()

######

## Density plot of counts for each sample

pdf(paste0(outGraphs, "/Normalize_and_Outliers_Density_preN_E10in80.pdf"))
plot(density(log2(fExDatDF[ ,1] + 1))
     , col = as.factor(metDatDF$Genotype.factor)[1]
     , main = paste(graphCodeTitle, "Density HTSeq log2(Counts + 1)"
                    , sep = "\n")
     , xlab = "Counts", ylim = c(0, 0.2))
for (i in 2:dim(fExDatDF)[2]) {
  lines(density(log2(fExDatDF[ ,i] + 1))
        , col = as.factor(metDatDF$Genotype.factor)[i])
}
dev.off()

######

## PC Analysis of Technical Covariates

# Select covariates
pairsDat <- cbind(metDatDF[c("RNAconcentration"
                             , "RIN"
                             , "BioanalyzerConcentration"
                             , "Mouse"
                             , "Genotype"
                             , "Sex"
                             , "DOB"
                             , "Cage"
                             , "Breeders"
                             , "DissectionDate"
                             , "DissectionTime"
                             , "DissectionOrder"
                             , "Brain"
                             , "ChipIndex"
                             , "GeneLengthScore"
                             , "GCcontent")]
                  , seqQCdatDF[c("TOTAL_READS"
                                 , "PF_READS_ALIGNED"
                                 , "MEDIAN_5PRIME_TO_3PRIME_BIAS"
                                 , "READ_PAIR_DUPLICATES"
                                 , "READ_PAIR_OPTICAL_DUPLICATES")])
# Designate factors
pairsDat[c("Genotype"
           , "Sex"
           , "DOB"
           , "Cage"
           , "Breeders"
           , "DissectionDate"
           , "DissectionTime"
           , "DissectionOrder"
           , "Brain")] <- lapply(pairsDat[c("Genotype"
                                            , "Sex"
                                            , "DOB"
                                            , "Cage"
                                            , "Breeders"
                                            , "DissectionDate"
                                            , "DissectionTime"
                                            , "DissectionOrder"
                                            , "Brain"
           )], as.factor)

pairsDat[c("TOTAL_READS"
           , "PF_READS_ALIGNED"
           , "READ_PAIR_DUPLICATES"
           , "READ_PAIR_OPTICAL_DUPLICATES")] <- lapply(pairsDat[
             c("TOTAL_READS"
               , "PF_READS_ALIGNED"
               , "READ_PAIR_DUPLICATES"
               , "READ_PAIR_OPTICAL_DUPLICATES"
             )], as.numeric)

head(seqQCdatDF)

# Run PCA

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(log(t(as.matrix(exDatDF) + 1), 2), scale = F))
topPCs <- Run_PCA(meanCenteredM)
Plot_Corr_Matrix(pairsDat, topPCs, title = paste(graphCodeTitle, "log2(reads + 1)", sep = "\n")
                 , outPath = paste(outGraphs, "/Normalize_and_Outliers_PCA_Covariates_E10in80_Log2readsP1.pdf", sep = ""))

######

# CNTNAP2 expression KO vs WT
cnWtDF <- data.frame(
  Count = t(fExDatDF["ENSMUSG00000039419"
                         , as.factor(metDatDF$Genotype.factor) == "WT"])[ ,1]
  , Genotype = "WT"
  , BrainRegion = metDatDF$Brain[as.factor(metDatDF$Genotype.factor) == "WT"])
cnKoDF <- data.frame(
  Count = t(fExDatDF["ENSMUSG00000039419"
                           , as.factor(metDatDF$Genotype.factor) == "KO"])[ ,1]
  , Genotype = "KO"
  , BrainRegion = metDatDF$Brain[as.factor(metDatDF$Genotype.factor) == "KO"])
cntnap2DF <- rbind(cnWtDF, cnKoDF)

ggDF <- melt(cntnap2DF)
ggplot(ggDF, aes(x = Genotype, y = log(value + 1, 2))) + 
  geom_boxplot(outlier.shape = NA) + #avoid plotting outliers twice
  geom_jitter(position = position_jitter(width = .4, height = 0), shape = 1, aes(col = BrainRegion)) +
  xlab("Genotype") +
  ylab("log2(Counts + 1)") +
  labs(title = paste(graphCodeTitle, "CNTNAP2 Expression: HTSeq log2(Counts + 1)"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(paste0(outGraphs, "/Normalize_and_Outliers_CNTNAP2_preN_E10in80.pdf"))
################################################################################

# (4) Normalize

# Use CQN to normalize for GC content and gene length and quantile
# Outputs values on log2 scale
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
                , sqn = FALSE)
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
ggsave(file = paste0(outGraphs, "/Normalize_and_Outliers_GC_Length_PostCQN_E10in80.pdf")
       , height = 6)
################################################################################

# (5) Remove Outliers

# We will go with the less rigorous 'All Samples' method here, but the 'By Dx' method is also included for your reference

## For All Samples

pdf("./SampleNetwork_allSamps.pdf")


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

metOrDatDF <- metDatDF[!flag, ]
cqnFexDat <- cqnFexDat[ ,!flag]
seqQcOrdatDF <- seqQCdatDF[!flag,]

# Format CQN normalized and outlier removed expression data so that metadata
# rows and expression data columns are in same sample order
cqnFexDat <- as.data.frame(cqnFexDat)
cqnFexDat <- cqnFexDat[ ,order(colnames(cqnFexDat))]
################################################################################

### (7) View Data Post CQN and Post Outlier Removal

## PC Analysis of Technical Covariates

# Select covariates
pairsDat <- cbind(metOrDatDF[c("RNAconcentration"
                             , "RIN"
                             , "BioanalyzerConcentration"
                             , "Mouse"
                             , "Genotype"
                             , "Sex"
                             , "Brain"
                             , "ChipIndex"
                             , "GeneLengthScore"
                             , "GCcontent")]
                  , seqQcOrdatDF[c("TOTAL_READS"
                                 , "PF_READS_ALIGNED"
                                 , "MEDIAN_5PRIME_TO_3PRIME_BIAS"
                                 , "READ_PAIR_DUPLICATES"
                                 , "READ_PAIR_OPTICAL_DUPLICATES")])
# Designate factors
pairsDat[c("Genotype"
           , "Sex"
           , "Brain")] <- lapply(pairsDat[c("Genotype"
                                            , "Sex"
                                            , "Brain"
           )], as.factor)

pairsDat[c("TOTAL_READS"
           , "PF_READS_ALIGNED"
           , "READ_PAIR_DUPLICATES"
           , "READ_PAIR_OPTICAL_DUPLICATES")] <- lapply(pairsDat[
             c("TOTAL_READS"
               , "PF_READS_ALIGNED"
               , "READ_PAIR_DUPLICATES"
               , "READ_PAIR_OPTICAL_DUPLICATES"
             )], as.numeric)

# Run PCA

# Do a PCA of the sequencing statistics

# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(as.matrix(cqnFexDat)), scale = F))
topPCsEx <- Run_PCA(meanCenteredM, "Expression")
meanCenteredM <- t(scale(seqQcOrdatDF[ ,-1], scale = F))
topPCsSeq <- Run_PCA(meanCenteredM, "SeqStat")
Plot_Corr_Matrix(cbind(pairsDat, topPCsSeq, topPCsEx), pairsDat$Brain
 , title = paste(graphCodeTitle, "CQN Normalized, Outliers Removed", sep = "\n")
 , outPath = paste(outGraphs
      , "/Normalize_and_Outliers_PCA_Covariates_E10in80_CQN_OR.pdf", sep = ""))

################################################################################

### (8) Regress out Sequencing Statistics PC2 and GC content, check differential
###     gene expression

meanCenteredM <- t(scale(t(as.matrix(cqnFexDat)), scale = F))
topPCsEx <- Run_PCA(meanCenteredM, "Expression")
meanCenteredM <- t(scale(seqQcOrdatDF[ ,-1], scale = F))
topPCsSeq <- Run_PCA(meanCenteredM, "SeqStat")

modelDat <- data.frame(Genotype = metOrDatDF$Genotype, SeqPC2 = topPCsSeq[ ,2]
                       , GCcontent = metOrDatDF$GCcontent)

design <- "~Genotype+SeqPC2+GCcontent"
mod <- model.matrix(as.formula(design), data = modelDat)
colnames(mod) <- make.names(colnames(mod))
lm <- lmFit(cqnFexDat, mod)
lm <- eBayes(lm)
fit <- lm

regCqnFexDat <- cqnFexDat
for (i in 1:nrow(cqnFexDat)){
  regCqnFexDat[i, ] <- cqnFexDat[i, ] - lm$coefficients[i,3] * mod[ ,3] - lm$coefficients[i,4] * mod[ ,4]
}

## Differential Gene Expression

# DGE for CQN, Outliers Removed, Regressed 
modelDat <- data.frame(Genotype = metOrDatDF$Genotype, SeqPC2 = topPCsSeq[ ,2]
                       , GCcontent = metOrDatDF$GCcontent)

design <- "~Genotype+SeqPC2+GCcontent"
mod <- model.matrix(as.formula(design), data = modelDat)
colnames(mod) <- make.names(colnames(mod))
lm <- lmFit(cqnFexDat, mod)
lm <- eBayes(lm)
fit <- lm
gt <- topTable(lm, coef = 2,number = nrow(regCqnFexDat), sort.by = "none")
hist(gt$adj.P.Val)
hist(gt$logFC, breaks = 100)
hist(gt$logFC, breaks = 100, xlim = c(-2, 2))

gt["ENSMUSG00000039419", ]
head(gt[order(gt$adj.P.Val), ], 50)
head(gt[order(-gt$logFC), ], 50)

# DGE for CQN, Outliers Removed
lm <- lmFit(cqnFexDat, mod)
lm <- eBayes(lm)
fit <- lm
gt <- topTable(lm, coef = 2,number = nrow(cqnFexDat), sort.by = "none")

head(gt[order(gt$adj.P.Val), ], 50)
head(gt[order(-gt$logFC), ], 50)
hist(gt$adj.P.Val)
hist(gt$logFC, breaks = 100)

# DGE for raw counts
modelDat <- data.frame(Genotype = metDatDF$Genotype)
design <- "~Genotype"
mod <- model.matrix(as.formula(design), data = modelDat)

lm <- lmFit(exDatDF, mod)
lm <- eBayes(lm)
fit <- lm
gt <- topTable(lm, coef = 2, number = nrow(exDatDF), sort.by = "none")

head(gt[order(gt$adj.P.Val), ], 50)
head(gt[order(-gt$logFC), ], 50)
hist(gt$adj.P.Val, breaks = 100)
hist(gt$logFC, breaks = 100)
hist(gt$logFC, breaks = 100, xlim = c(-200, 200))
################################################################################

### (9) View Data Post CQN and Post Outlier Removal

## PC Analysis of Technical Covariates

# PCA of Expression Values
# Centers the mean of all genes - this means the PCA gives us the eigenvectors
# of the geneXgene covariance matrix, allowing us to assess the proportion of
# variance each component contributes to the data
meanCenteredM <- t(scale(t(as.matrix(regCqnFexDat)), scale = F))
topPCsEx <- Run_PCA(meanCenteredM, "Expression")

# PCA of the sequencing statistics
meanCenteredM <- t(scale(seqQcOrdatDF[ ,-1], scale = F))
topPCsSeq <- Run_PCA(meanCenteredM, "SeqStat")

Plot_Corr_Matrix(cbind(pairsDat, topPCsSeq, topPCsEx), pairsDat$Brain
 , title = paste(graphCodeTitle, "CQN Normalized, Outliers Removed, SeqStats PC2 and 3 Regressed out", sep = "\n")
 , outPath = paste(outGraphs
  , "/Normalize_and_Outliers_PCA_Covariates_E10in80_CQN_OR_Rg.pdf", sep = ""))

## MDS All Genes

# Calculate MDS
mdsDF <- calcMDS(regCqnFexDat)
# Add column with sample type info
mdsDF$type <- metOrDatDF$Genotype.factor

ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 1) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1, size = 1.5) +
  scale_color_discrete(name = "Sample Type") +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq + CQN + Outliers Removed + Reg"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(paste0(outGraphs, "/Normalize_and_Outliers_MDS_GT_E10in80_CQN_OR_Rg.pdf"))

mdsDF$type <- metOrDatDF$Brain
ggplot(mdsDF, aes(x = X1, y = X2, color = factor(type))) +
  geom_point(size = 1) +
  geom_text(aes(label = row.names(mdsDF)), vjust = -1, size = 1.5) +
  scale_color_discrete(name = "Brain Region") +
  xlab(paste("PC1 (", signif(100*mdsDF$pc1, 3), "%)", sep = "")) +
  ylab(paste("PC2 (", signif(100*mdsDF$pc2, 3), "%)", sep = "")) +
  labs(title = paste(graphCodeTitle, "MDS Plot: HTSeq + CQN + Outliers Removed + Reg"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black")) +
  theme(aspect.ratio = 4/5)
ggsave(paste0(outGraphs, "/Normalize_and_Outliers_MDS_BR_E10in80_CQN_OR_Rg.pdf"))

#####

## Boxplot of counts for each sample
pdf(paste0(outGraphs, "/Normalize_and_Outliers_BoxPlot_E10in80_CQN_OR_Rg.pdf"))
par(las = 2)
par(mar = c(14,4,3,1))
boxplot(regCqnFexDat, range = 0, col = as.factor(metOrDatDF$Genotype.factor)
        , main = paste(graphCodeTitle, "Boxplot HTSeq + CQN + Outliers Removed + Reg"
                       , sep = "\n")
        , xlab = "", ylab = "Counts")
par(las = 1)
mtext("Sample", side = 1, line = 12)
dev.off()

#####

## Density plot of counts for each sample
pdf(paste0(outGraphs, "/Normalize_and_Outliers_Density_E10in80_CQN_OR_Rg.pdf"))
plot(density(regCqnFexDat[ ,1])
     , col = as.factor(metOrDatDF$Genotype.factor)[1]
     , main = paste(graphCodeTitle, "Density HTSeq + CQN + Outliers Removed + Reg"
                    , sep = "\n")
     , xlab = "Counts", ylim = c(0, 0.01))
for (i in 2:dim(cqnFexDat)[2]) {
  lines(density(regCqnFexDat[ ,i])
        , col = as.factor(metOrDatDF$Genotype.factor)[i])
}
dev.off()

#####

## Boxplot CNTNAP2 expression KO vs WT

cnWtDF <- data.frame(
  Count = t(cqnFexDat["ENSMUSG00000039419"
                     , as.factor(metOrDatDF$Genotype.factor) == "WT"])[ ,1]
  , Genotype = "WT"
  , BrainRegion = metOrDatDF$Brain[as.factor(metOrDatDF$Genotype.factor) == "WT"])
cnKoDF <- data.frame(
  Count = t(cqnFexDat["ENSMUSG00000039419"
                     , as.factor(metOrDatDF$Genotype.factor) == "KO"])[ ,1]
  , Genotype = "KO"
  , BrainRegion = metOrDatDF$Brain[as.factor(metOrDatDF$Genotype.factor) == "KO"])
cntnap2DF <- rbind(cnWtDF, cnKoDF)

ggDF <- melt(cntnap2DF)
ggplot(ggDF, aes(x = Genotype, y = log(value + 1, 2))) + 
  geom_boxplot(outlier.shape = NA) + #avoid plotting outliers twice
  geom_jitter(position = position_jitter(width = .4, height = 0), shape = 1, aes(col = BrainRegion)) +
  xlab("Genotype") +
  ylab("Expression: CQN Normalized") +
  labs(title = paste(graphCodeTitle, "CNTNAP2 Expression: HTSeq + CQN + Outliers Removed + Reg"
                     , sep = "\n")) +
  theme_grey(base_size = 16) +
  theme(axis.text = element_text(color = "black"))
ggsave(paste0(outGraphs, "/Normalize_and_Outliers_CNTNAP2_E10in80_CQN_OR_Rg.pdf"))

save(cqnFexDat, metOrDatDF, file = "../data/Expr_E10in80_Nm_OR_Rg.rda")
################################################################################
