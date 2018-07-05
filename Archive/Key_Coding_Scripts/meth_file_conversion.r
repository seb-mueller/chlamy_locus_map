#Code written by Nick Matthews to make txt files of methylation data from Tom Hardcastle into gff3 files for use in annotation
#Data: 19/01/16


setwd("/home/bioinf/nem34/segmentMap_II/meth_data/")

library(rtracklayer)

CG <- read.delim("chlamy_CGmeth.txt",header=TRUE)
CG2 <- GRanges(seqnames=CG$seqnames,ranges = IRanges(start=CG$start,end=CG$end),strand=CG$strand)
CG2$ordering <- rep("none", length(CG2))
export.gff3(CG2,"chlamy_CGmeth.gff3")
#import.gff3("chlamy_CGmeth.gff3")

CHH <- read.table("chlamy_CHHmeth.txt",header=TRUE)
CHH2 <- GRanges(seqnames=CHH$seqnames,ranges = IRanges(start=CHH$start,end=CHH$end),strand=CHH$strand)
CHH2$ordering <- rep("none", length(CHH2))
export.gff3(CHH2,"chlamy_CHHmeth.gff3")
#import.gff3("chlamy_CHHmeth.gff3")

CHG <- read.table("chlamy_CHGmeth.txt",header=TRUE)
CHG2 <- GRanges(seqnames=CHG$seqnames,ranges = IRanges(start=CHG$start,end=CHG$end),strand=CHG$strand)
CHG2$ordering <- rep("none", length(CHG2))
export.gff3(CHG2,"chlamy_CHGmeth.gff3")
#import.gff3("chlamy_CHGmeth.gff3")

diffCG <- read.table("Andy_differentialMeth_CG.txt",header=TRUE)
diffCG2 <- GRanges(seqnames=diffCG$seqnames,ranges = IRanges(start=diffCG$start,end=diffCG$end),strand=diffCG$strand,ordering=diffCG$ordering)
export.gff3(diffCG2,"Andy_differentialMeth_CG.gff3")
#import.gff3("Andy_differentialMeth_CG.gff3")

diffCHH <- read.table("Andy_differentialMeth_CHH.txt",header=TRUE)
diffCHH2 <- GRanges(seqnames=diffCHH$seqnames,ranges = IRanges(start=diffCHH$start,end=diffCHH$end),strand=diffCHH$strand,ordering=diffCHH$ordering)
export.gff3(diffCHH2,"Andy_differentialMeth_CHH.gff3")
#x<-import.gff3("Andy_differentialMeth_CHH.gff3")

diffCHG <- read.table("Andy_differentialMeth_CHG.txt",header=TRUE)
diffCHG2 <- GRanges(seqnames=diffCHG$seqnames,ranges = IRanges(start=diffCHG$start,end=diffCHG$end),strand=diffCHG$strand,ordering=diffCHG$ordering)
export.gff3(diffCHG2,"Andy_differentialMeth_CHG.gff3")
#import.gff3("Andy_differentialMeth_CHG.gff3")
