## Processing the QC files
## Obtain .csv's summed across samples with 4.5_CompileQC.sh

setwd("~Path/to/your/directory/QC")
options(stringsAsFactors=FALSE)

pdf("./PicardQC.pdf")
picard.align <- read.table("./RNAseqAlign.txt",header=TRUE)
rownames(picard.align) <- picard.align[,1]
picard.align <- picard.align[,-c(1)]

colnames(picard.align) <- c("CATEGORY","TOTAL_READS","PF_READS","PCT_PF_READS","PF_NOISE_READS",
                            "PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PF_ALIGNED_BASES","PF_HQ_ALIGNED_READS",
                            "PF_HQ_ALIGNED_BASES","PF_ HQ_ALIGNED_Q20_BASES","PF_HQ_MEDIAN_MISMATCHES",
                            "PF_MISMATCH_RATE","PF_HQ_ERROR_RATE","PF_INDEL_RATE","MEAN_READ_LENGTH",
                            "READS_ALIGNED_IN_PAIRS","PCT_READS_ALIGNED_IN_PAIRS","BAD_CYCLES","STRAND_BALANCE",
                            "PCT_CHIMERAS","PCT_ADAPTER")

summary(picard.align)

perc.hq <- picard.align[,"PF_HQ_ALIGNED_BASES"]/picard.align[,"PF_ALIGNED_BASES"]
quantile(perc.hq)

## Transcript level QC
picard.rnaseq.qc <- read.table("./RNAseqQC.txt",header=TRUE)
#colnames(picard.rnaseq.qc) <- c("FileName","PF_BASES","PF_ALIGNED_BASES","CODING_BASES","UTR_BASES","INTRONIC_BASES","INTERGENIC_BASES","IGNORED_READS","CORRECT_STRAND_READS","INCORRECT_STRAND_READS",
#                                                                "PCT_CODING_BASES","PCT_UTR_BASES","PCT_INTRONIC_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES","PCT_USABLE_BASES",
#                                                                "PCT_CORRECT_STRAND_READS","MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS","MEDIAN_5PRIME_TO_3PRIME_BIAS")

summary(picard.rnaseq.qc)

quantile(picard.rnaseq.qc[,"MEDIAN_5PRIME_TO_3PRIME_BIAS"],c(seq(0,1,0.05)))
quantile(picard.rnaseq.qc[,"MEDIAN_5PRIME_TO_3PRIME_BIAS"],c(0.025,0.975))

quantile(picard.rnaseq.qc[,"PF_BASES"]/100) ## gives the number of fragments aligned - 38.9 million is the median
quantile(picard.rnaseq.qc[,"CODING_BASES"]/100) ## gives the number of fragments aligned to coding regions - 8.5 million on average
quantile(picard.rnaseq.qc[,"INTRONIC_BASES"]/100) ## gives the number of fragments aligned to intronic regions - 15 million on average
quantile(picard.rnaseq.qc[,"INTERGENIC_BASES"]/100) ## gives the number of fragments aligned to intronic regions - 2.7 million on average   



tx.coverage <- read.table("./TranscriptCoverage.txt") 
keep <- match(unique(tx.coverage[,1]),tx.coverage[,1])
tx.coverage <- tx.coverage[keep,]
rownames(tx.coverage) <- tx.coverage[,1]
tx.coverage.cols <- as.character(tx.coverage[1,seq(4,204,by=2)])
tx.coverage <- tx.coverage[,seq(5,205,by=2)]
colnames(tx.coverage) <- tx.coverage.cols
tx.coverage.quant <- apply(tx.coverage,2,quantile,c(0.025,0.5,0.975))

par(mfrow=c(2,1))
boxplot(cbind(perc.hq,picard.rnaseq.qc[,c(15,11:14)]),ylab="Percent of bases",main="RNA-seq QC metrics across samples",names=c("High\nQuality","mRNA\nBases","Protein\nCoding","Untranslated\nRegion","Intronic\nRegion","Intergenic\nRegion"),ylim=c(0,1))

plot(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[2,],xlab="Percentile of gene body (5' -> 3')",ylab="Coverage relative to whole transcript",pch=19,ylim=c(0,1.7),main="Relative transcript coverage - median with 95% CIs across samples")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[2,],col="black")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[1,],col="grey")
lines(x=as.numeric(colnames(tx.coverage.quant)),y=tx.coverage.quant[3,],col="grey")


## GC content
## Columns c("GC","WINDOWS","READ_STARTS","MEAN_BASE_QUALITY","NORMALIZED_COVERAGE","ERROR_BAR_WIDTH")
tmp <- read.table("./RNAseqGC.txt",skip=1)
picard.rnaseq.gc <- matrix(NA,nrow=nrow(tmp),ncol=101)
colnames(picard.rnaseq.gc) <- as.character(tmp[1,seq(2,602,by=6)])
rownames(picard.rnaseq.gc) <- tmp[,1]
picard.rnaseq.gc <- tmp[,seq(6,607,by=6)]
gc.coverage.quant <- apply(picard.rnaseq.gc,2,quantile,c(0.025,0.5,0.975))
colnames(gc.coverage.quant) <- as.character(tmp[1,seq(2,602,by=6)])

picard.rnaseq.gcqual <- matrix(NA,nrow=nrow(tmp),ncol=101)
colnames(picard.rnaseq.gcqual) <- tmp[1,seq(2,602,by=6)]
rownames(picard.rnaseq.gcqual) <- tmp[,1]
picard.rnaseq.gcqual <- tmp[,seq(5,607,by=6)]
gc.qual.quant <- apply(picard.rnaseq.gcqual,2,quantile,c(0.025,0.5,0.975))

picard.rnaseq.gcbin <- matrix(NA,nrow=nrow(tmp),ncol=101)
colnames(picard.rnaseq.gcbin) <- tmp[1,seq(2,602,by=6)]
rownames(picard.rnaseq.gcbin) <- tmp[,1]
picard.rnaseq.gcbin <- tmp[,seq(3,607,by=6)]
gc.bins.quant <- apply(picard.rnaseq.gcbin,2,quantile,c(0.5))
gc.bins.quant <- gc.bins.quant/sum(gc.bins.quant)

par(mfrow=c(3,1))
plot(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[2,],xlab="% GC content",ylab="Normalized coverage",pch=19,ylim=c(0,15),main="Coverage by GC percentage - median with 95% CIs across samples")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[2,],col="black")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[1,],col="grey")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.coverage.quant[3,],col="grey")

plot(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[2,],xlab="%GC content",ylab="Read quality score",pch=19,ylim=c(0,36),main="Read quality by GC percent - median with 95% CIs across samples")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[2,],col="black")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[1,],col="grey")
lines(x=as.numeric(colnames(gc.coverage.quant)),y=gc.qual.quant[3,],col="grey")

plot(density(x=as.numeric(colnames(gc.coverage.quant)),y=gc.bins.quant),xlab="%GC content",ylab="Proportion of 100bp bins",pch=19,ylim=c(0,max(gc.bins.quant)+0.01),main="Proportion of bins corresponding to each GC percentile")

##Duplication Metrics
## Alter RNAseqDuplication.txt to a .csv manually first

picard.rnadups.qc <- read.csv(file="./RNAseqDuplication.csv")
summary(picard.rnadups.qc)

## More GC metrics
gc.summary <- read.table("./GCsummary.txt")
colnames(gc.summary)=as.character(gc.summary[1,])
gc.summary=gc.summary[-1,]
rownames(gc.summary)=as.character(gc.summary[,1])
gc.summary=gc.summary[,-1]
gc.summary = sapply(gc.summary[,1:5],as.factor,USE.NAMES = FALSE)
summary(gc.summary)

## Compile most important QC metrics - choose the metrics that you think are valuable to the analysis and to sequencing statistics
## Metrics concering read depth (ie. Total Reads) and GC bias (ie. 5'/3' bias) should always be included

PQCdat <- cbind(picard.align[,c(2,6,9,21)],picard.rnaseq.qc,
                                gc.summary[,4],picard.rnadups.qc[,c(7,8)])
write.csv(PQCdat,"./PicardToolsQC.csv")
dev.off()
