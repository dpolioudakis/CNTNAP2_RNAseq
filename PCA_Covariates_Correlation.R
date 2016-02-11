# Damon Polioudakis
# 2016-02-09
# Compare technical covariates to principal components of gene expression data

# Inputs
#   HTseq counts expression
#   Metadata
#   Picard QC
#   GC bias and Gene length bias

# Outputs
#   Correlation plot of principal compenents and covariates
################################################################################

rm(list=ls())
sessionInfo()

library(xlsx)
library(boot)
library(reshape2)

# Load data and assign variables

# Input data
# Gene expression from HTSeq
exDatDF <- read.csv("../data/HTSC/Exprs_HTSCexon.csv", row.names = 1)
# Metadata
metDatDF <- read.xlsx("../metadata/Cntnap2 WT KO RNA-seq analysis.xlsx"
                      , sheetIndex = 1)
# Picard QC data
seqQCdatDF <- read.csv("../metadata/PicardToolsQC.csv")
# Gene Length and GC bias
load("../analysis/tables/Gene_Length_GC_Bias.rda")
lenBiasDF <- data.frame(GeneLengthScore = lenBias)
gcBiasDF <- data.frame(GCcontent = gcBias)

# Output file path prefix
outPathCorr <- "../analysis/graphs/PCA_Covariates_Correlation"
################################################################################

# Format data sets

# Add length bias to metadata
row.names(lenBiasDF) <- substring(row.names(lenBiasDF), 2)
metDatDF <- merge(x = metDatDF, y = lenBiasDF
                  , by.x = "Sample.ID", by.y = "row.names")
# Add GC bias to metadata
row.names(gcBiasDF) <- substring(row.names(gcBiasDF), 2)
metDatDF <- merge(x = metDatDF, y = gcBiasDF
                  , by.x = "Sample.ID", by.y = "row.names")

# Filter for genes exp > 10 counts in  more than 80% of samples
pres <- apply(exDatDF > 10, 1, sum)
idx <- which(pres > 0.8 * dim(exDatDF)[2])
fExDatDF <- exDatDF[idx, ]
################################################################################

# PC Analysis

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
Plot_Corr_Matrix <- function (pairsDat, topPCs, title, outPath) {
  pdf(outPath, height = 20, width = 24)
  pairs(cbind(pairsDat, topPCs), pch = 19
        , upper.panel = panel.cor
        , main = paste("Covariates and MaxQuant Comparison -- |Spearman's rho| correlation values",
                       "\n", title, sep = "")
        # , cex.labels = (cex = 1.25)
  )
  dev.off()
}

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
meanCenteredM <- t(scale(log(t(as.matrix(exDatDF[1:(nrow(exDatDF)-97), ])+1), 2)
                         , scale=F))
topPCs <- Run_PCA(meanCenteredM)
Plot_Corr_Matrix(pairsDat, topPCs, title = "log2(reads + 1)"
                 , outPath = paste(outPathCorr, "_Log2readsP1.pdf", sep = ""))







