# Damon Polioudakis
# 2016-02-08
# Calculate gene length and GC content for exon union of each gene in list
################################################################################


##### INCORRECT OUTPUT - NEEDS DEBUGGING AT INDICATED AREA ######
# Check that union exon lengths and gc correlate with longest isoform lengths
# and gc

rm(list=ls())
sessionInfo()

source("http://bioconductor.org/biocLite.R")
library(Repitools)
library(BSgenome)
# Run if BSgenome.Mmusculus.UCSC.mm10 is not installed:
# biocLite("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome.Mmusculus.UCSC.mm10)
library(GenomicRanges)
# Run if Genominator is not installed:
# biocLite("Genominator")
library(Genominator)

options(stringsAsFactors=FALSE)

# Load data and assign variables

exDatDF <- read.csv("../data/HTSC/Exprs_HTSCexon.csv")

# Get Gencode 18 gtf file - this was cleaned by selecting only the columns
# containing the word "exon" and with the relevant information - chrnum,
# feature type, strand, start, end, and ensg and ense IDs separated by a semicolon
gtfInfoDF <- read.table("../mouseRefGenome/gencode.vM8.annotation.gtf", sep = "\t")
################################################################################

# Format and Filter

# Keep only the exon level features
keep <- gtfInfoDF[ , 3]=="exon"
gtfInfoDF <- gtfInfoDF[keep, ]

# Split the semicolon separated information
geneExonInfo <- unlist(strsplit(gtfInfoDF[ , 9], "[;]"))

# Finding which has gene_id for ensembl ID
genCol<- which(regexpr("gene_id ", geneExonInfo) > 0)
getSeq <- geneExonInfo[genCol]
ENSGID <- substr(getSeq, 9, 100)
length(unique(ENSGID)) #46983

transCol <- which(regexpr("transcript_id ", geneExonInfo) > 0)
tranSeq <- geneExonInfo[transCol]
ENSEID <- substr(tranSeq, 16, 100)
length(unique(ENSEID)) #113231
gtfInfoDF <- cbind(gtfInfoDF[ , c(1:8)], ENSGID, ENSEID)

geneCol <- which(regexpr("gene_name ", geneExonInfo) > 0)
geneSeq <- geneExonInfo[geneCol]
GENEID <- substr(geneSeq, 12, 100)
length(unique(GENEID)) #46882

gtfInfoDF <- cbind(gtfInfoDF[ , c(1:8)], ENSGID, ENSEID)
# 6 and 8 columns are blank - remove
gtfInfoDF <- gtfInfoDF[ , -c(6, 8)]


################ FIX - remove? ############ -> this is only keeping 1st exon?
# Keep only one copy of each ENSEID - the gtf file records one copy for each
# transcript id
keep <- match(unique(ENSEID), ENSEID)
gtfInfoDF1 <- gtfInfoDF[keep, ]
##gtfInfoDF[,1] <- substr(gtfInfoDF[,1],4,10) ## 672406 exons is exactly what biomaRt contains
##########################################


# Recode things for the Genominator package
# Using as.factor to coerce chromosome names can really botch things up.. beware!
# So go ahead and convert MT, X, and Y to numbers throughout, unless necessary
# for other purposes
chrNums <- gtfInfoDF1[ , 1]
chrNums[chrNums=="chrM"] <- "20"
chrNums[chrNums=="chrX"] <- "21"
chrNums[chrNums=="chrY"] <- "22"
# Convert UCSC Chromsome names to Ensembl (chr1 -> 1)
chrNums <- gsub("chr", "", chrNums)
# rmChR.col1=which(regexpr("HG", chrNums)>0)
# rmChR.col2= which(regexpr("GL", chrNums)>0) ## removing Non-annotated(NT) chromosomes
# rmChR.col3= which(regexpr("HS", chrNums)>0)
# rmChR.col=c(rmChR.col1,rmChR.col2,rmChR.col3)
# gtfInfoDF1=gtfInfoDF1[-rmChR.col,]
# chrNums <- chrNums[-rmChR.col]
gtfInfoDF1[ ,1] <- chrNums ## Check here

gtfInfoDF <- gtfInfoDF1

strinfo <- as.character(gtfInfoDF[ , 6])
strinfo[strinfo=="+"] <- 1L
strinfo[strinfo=="-"] <- -1L
gtfInfoDF[ , 6] <- strinfo

# chr integer, strand integer (-1L,0L,1L), start integer, end integer, ensg
# and transcript id
geneDat1 <- gtfInfoDF[ , c(1, 6, 4, 5, 7, 8)]
geneDat1 <- data.frame(as.numeric(chrNums), as.numeric(geneDat1[ , 2])
                       , as.numeric(geneDat1[ , 3]) ,as.numeric(geneDat1[ , 4])
                       , geneDat1[ , 5] ,geneDat1[ , 6])
names(geneDat1) <- c("chr","strand","start","end","ensembl_gene_id","ensembl_exon_id")

geneDatX <- geneDat1[order(geneDat1[ , 1],geneDat1[ , 3]), ]
# Remove NAs from ERCC chromosomes
# geneDatX <- geneDatX[complete.cases(geneDatX),]
# Have genominator check if this is a valid data object
validAnnotation(geneDatX)
# Should take few minutes !!!
geneDatX <- makeGeneRepresentation(annoData = geneDat1, type = "Ugene"
                                   , gene.id = "ensembl_gene_id"
                                   , transcript.id = "ensembl_exon_id"
                                   , verbose = TRUE) 

save(geneDatX, file = "../source/Genominator_Union_GeneModels_m38.rda")
load(file = "../source/Genominator_Union_GeneModels_m38.rda")

# Now use the genominator output to calculate GC content ### Use mac laptop ###
geneDat2 <- cbind(geneDatX, geneDatX[ , 3] - geneDatX[ , 2])
geneDat2 <- geneDat2[order(geneDat2[ , 5]) , ]

# Change formatting again
chrNums <- geneDat2[ ,"chr"]
chrNums[chrNums=="20"] <- "M" ## important as UCSC codes MT as M
chrNums[chrNums=="21"] <- "X"
chrNums[chrNums=="22"] <- "Y"
stInfo <- geneDat2[ ,"strand"]
stInfo[stInfo==c(-1)] <- "-"
stInfo[stInfo==c(1)] <- "+"

# Calculate GC content from mm10 using the union exon ranges determined by
# selecting "exons" from the gtf file above
# Convert to a genomic ranges object
gcQuery <- GRanges(paste("chr", chrNums, sep = "")
                   , IRanges(geneDat2[ , 2], geneDat2[ , 3]), strand = stInfo)
gcContent <- gcContentCalc(x = gcQuery, organism = Mmusculus)

# Take a length weighted average of GC content percentages to get the GC content
# for the union gene model
head(geneDat2)
geneDat2 <- cbind(geneDat2, gcContent)
geneDat2 <- cbind(geneDat2, gcContent * geneDat2[ , 6])
unionGenes <- by(geneDat2[ , 6], as.factor(geneDat2[ , 5]), sum)
unionGC <- by(geneDat2[ , 8], as.factor(geneDat2[ , 5]), sum)
geneDat3 <- cbind(unionGenes, unionGC / unionGenes)
colnames(geneDat3) <- c("UnionGeneLength", "UnionGCcontent")
mm10exUnionAnno <- geneDat3

## Save for further usage
save(mm10exUnionAnno, file = "../source/Ensembl_Mm10_Exon_Union_Anno.rda")
load("../source/Ensembl_Mm10_Exon_Union_Anno.rda")

# Calculate exon length bias
unionGenes <- data.frame(Length = mm10exUnionAnno[ , 1])
exLenDF <- merge(x = exDatDF, y = unionGenes, by.x = "X", by.y = "row.names" )
lenBias <- apply(exLenDF, 2
                 , function(counts) sum(as.numeric(counts) * exLenDF["Length"]) / 
                   sum(as.numeric(counts)))
gcBiasDF <- merge(x = exDatDF, y = mm10exUnionAnno
                  , by.x = "X", by.y = "row.names" )
gcBiasDF <- gcBiasDF[complete.cases(gcBiasDF), ]
# Select only columns with counts to calculate GC content bias
gcBias <- apply(gcBiasDF[ ,-c(1, (ncol(gcBiasDF)-1):ncol(gcBiasDF))], 2
                , function(counts) sum(as.numeric(counts)
                                       * gcBiasDF["UnionGCcontent"])
                                       / sum(as.numeric(counts)))
lenBias <- lenBias[-c(1, length(lenBias))]
save(lenBias, gcBias, file = "../analysis/tables/Gene_Length_GC_Bias.rda")







## Check difference in relationship to GC before and after normalization
preNorm <- datExpr.HTSC
postNorm <- cqn.dat
keepgenes <- intersect(rownames(preNorm),rownames(postNorm))
preNorm <- preNorm[match(keepgenes,rownames(preNorm)),]
postNorm <- postNorm[match(keepgenes,rownames(postNorm)),]
geneAnno1 <- geneAnno[match(keepgenes,rownames(geneAnno)),]

qualMat <- matrix(NA,nrow=ncol(preNorm),ncol=4)
colnames(qualMat) <- c("pre.Norm.GC.cor","pre.Norm.Length.cor","post.Norm.GC.cor","post.Norm.Length.cor")

for (i in 1:nrow(qualMat)) {
  qualMat[i,1] <- cor(preNorm[,i],geneAnno1[,2],method="spearman")
  qualMat[i,2] <- cor(preNorm[,i],geneAnno1[,1],method="spearman")
  qualMat[i,3] <- cor(postNorm[,i],geneAnno1[,2],method="spearman")
  qualMat[i,4] <- cor(postNorm[,i],geneAnno1[,1],method="spearman")
}
quantile(qualMat[,1])
quantile(qualMat[,3])
quantile(qualMat[,2])
quantile(qualMat[,4])
pdf("GCcorrelations.pdf",width=8,height=8)
par(mfrow=c(2,2))
hist(qualMat[,1],main=colnames(qualMat)[1],xlim=c(-0.15,0.1),ylim=c(0,150),breaks=seq(-0.15,0.1,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.15)
hist(qualMat[,3],main=colnames(qualMat)[3],xlim=c(-0.15,0.1),ylim=c(0,150),breaks=seq(-0.15,0.1,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.15)
hist(qualMat[,2],main=colnames(qualMat)[2],xlim=c(0,0.35),ylim=c(0,150),breaks=seq(0,0.35,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.10)
hist(qualMat[,4],main=colnames(qualMat)[4],xlim=c(0,0.35),ylim=c(0,150),breaks=seq(0,0.35,by=0.01),xlab="Spearman's rho across samples")
abline(v=0.10)
dev.off()