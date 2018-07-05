#Script to read in gene gff3 file and calculate intronic sequence
#Author: Nick Matthews
#Date: 11/06/18

#####Set-up#####
#Set working directory:
#Extract date
date <- format(Sys.Date(), "%d%m%Y")

#Load in necessary libraries
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(GenomicRanges)
library(rtracklayer)

#set working directory
setwd("C:/Users/Nick/Documents/Uni Work/Third Year/Project/segmentMap_II/Old_Annotation_Files")

#####Load in dataset#####
allGenes <- import.gff3("Creinhardtii_281_v5.5.gene_exons.gff3")

#Let's explore the data
widths <- data.frame(type = c(), width = c())
widths <- data.frame(type=c(as.character(unique(allGenes$type))),
                     width=sapply(unique(allGenes$type),function(x) mean(width(allGenes[allGenes$type==x]))))

#####Extract exonic sequences#####
genes <- allGenes[allGenes$type=="gene",]
exons <- allGenes[allGenes$type=="exon",]
mRNA <- allGenes[allGenes$type=="mRNA",]

#####Derive intronic sequences#####
#remove mRNA which don't have at least two exons overlapping them - therefore wouldn't have an intron by definition
exons <- unique(exons)
exonOverlaps <- data.frame(table(queryHits(findOverlaps(mRNA,exons))))
exonOverlaps$Keep <- rep(FALSE,nrow(exonOverlaps))
exonOverlaps$Keep[exonOverlaps$Freq > 1] <- TRUE
#apply to every gene?
mRNAList <- split(mRNA[exonOverlaps$Keep],1:length(mRNA[exonOverlaps$Keep]))
#Use endoapply to derive introns from areas of mRNA not covered by exons
intronList <- endoapply(mRNAList, function(x){
  introns <- NULL
  introns <- GenomicRanges::setdiff(x,exons)
  if(length(introns) > 0) {
    #If intron found, add in relevant information
    introns$source <- paste0("CalculatedFrom:",x$source)
    introns$parent <- x$ID
    introns$type <- "intron"
    introns$id <- paste0(x$ID,".intron.",1:length(introns))
    introns$pacid <- x$pacid
  } else {
    #If no introns, just return empty Granges object
    introns = GRanges()
  }
  introns
})

#Unlist into one GRanges object
introns <- unlist(intronList)

#####Save files#####
#Save introns as Rdata file and gff3 file
export.gff3(introns,paste0("CalculatedIntrons",date,".gff3"))

