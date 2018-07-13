#Script for running annotation functions
#Adapted by Nick Matthews from code written by Tom Hardcastle and Sebastian Muller
#Date: 09/02/16

#set working directory to github repository on cluster
setwd("/projects/nick_matthews/chlamy_locus_map_github")
#setwd("C:/Users/Nick/Documents/Uni Work/Third Year/Project/chlamy_locus_map")

#Load necesary libraries
library(xtable)
library(rtracklayer)
library(reshape)
library(segmentSeq)
library(pROC)
library(MASS)
library(RColorBrewer)
library(mclust)
library(baySeq)
library(MKmisc)
library(simpleboot)
source("./scripts/chlamy_source_code.r")

#Specify location of segmentation
segLocation <- "/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066"
#segLocation <- "C:/Users/Nick/Documents/Uni Work/Third Year/Project/segmentation_with_externals.r_2015-11-12_17/"
annoDir <- "/projects/nick_matthews/resources"

#####Load in and process segmentation, selecting significant loci#####
load(file.path(segLocation,"segD_first_chlamy_segmentation_nick.RData"))

# Select loci based on some fdr
loci7 <- selectLoci(segD, FDR = 0.1, perReplicate = TRUE) #

## creating/exporting coordinates object
gr7 <- loci7@coordinates 

# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr7) <- sprintf("CRSL%05.f0",1:length(gr7))

#load("gr7_just_cb.RData") #load gr7 you want to add to it...

# export as gff3 file for viewing in browser
export.gff3(gr7,con="loci7_fdr01.gff3")
#import.gff3("loci7_fdr01.gff3")

#####These next few functions compute and compile annotation files, this shouldn't need runnding every time#####
#Compute introns - this takes a while, don't run unless necessary
#intronCalculate()
#Process transposon file
transposonProcess(annoDir)
#Compute and compiles annotations
compileAnnotations(annoDir)
#####the next set of functions take the annotated locus object 'gr7' and add some more annotation data#####
#The functions are found 'chlamy_source_code.R' 

# annotate by size class
gr7 <- sizeClass(gr7,annoDir)

# annotate with overlapping features
gr7 <- featureAnn(gr7)

gr7 <- expressionClass(gr7)

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
cl <- makeCluster(24)
gr7 <- countingBiases(gr7,cl,segLocation)
stopCluster(cl)

#Old methylation function - may still be usefull, picks up a lot more methylation
gr7<-methylation1(gr7,annoDir)

#New methylation functions
gr7<-methylation2(gr7,annoDir)
gr7 <- methylationDiff(gr7annoDir) #Almost no results - very small datasets

#Extra Current annotations
gr7<-strainSpec(gr7,loci7)
gr7<-lifeCycle(gr7,loci7)
gr7<-mutantSpec(gr7,loci7)
gr7 <- phaseMatch(gr7)

#Other functions which were done for arabidopsis

#gr7 <- histoneAnnotate(gr7)
#cl <- makeCluster(24)
#gr7 <- methAnnotate(gr7, cl)
#stopCluster(cl)
#gr7 <- tissueSpec(gr7, loci7) #adapted for zygotes, key mutants, strains
#gr7 <- agoIP(gr7, loci7)
#cl <- makeCluster(24)
#Pol45(gr7, cl)
#stopCluster(cl)
#gr7 <- annPol(gr7)
#Save file
save(gr7,loci7,file="gr7_all.rdata")

