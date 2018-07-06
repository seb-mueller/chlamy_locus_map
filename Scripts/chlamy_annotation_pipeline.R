#Script for running annotation functions
#Adapted by Nick Matthews from code written by Tom Hardcastle and Sebastian Muller
#Date: 09/02/16

###### setting variables

fdr <- 0.05

# Specify location of segmentation
segLocation <- "/projects/nick_matthews/segmentation_2018/"
annoDir     <- "/projects/nick_matthews/resources"
#set working directory to github repository on cluster
gitdir      <- "/projects/nick_matthews/chlamy_locus_map_github"
setwd(gitdir)

cl <- makeCluster(24)
# using instead of arbitray versions
# e.g. "1f6085a"
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)

##### Load necesary libraries
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
source(file.path(gitdir, "/Scripts/chlamy_source_code.R"))

#####Load in and process segmentation, selecting significant loci#####
load(file.path(segLocation,"segD_first_chlamy_segmentation_nick.RData"))

# Select loci based on some fdr
loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE) #

## creating/exporting coordinates object
gr <- loci@coordinates

# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr) <- sprintf("CRSL%05.f0", 1:length(gr))

#load("gr_just_cb.RData") #load gr you want to add to it...

# export as gff3 file for viewing in browser
export.gff3(gr, con = file.path(segLocation, paste0("loci_fdr",myfdr, ".gff")))
#import.gff3("loci_fdr01.gff3")

#####These next few functions compute and compile annotation files, this shouldn't need runnding every time#####
#Compute introns - this takes a while, don't run unless necessary
#intronCalculate()
#Process transposon file
transposonProcess(annoDir)
#Compute and compiles annotations
compileAnnotations(annoDir)
#####the next set of functions take the annotated locus object 'gr' and add some more annotation data#####
#The functions are found 'chlamy_source_code.R'

# annotate by size class
gr <- sizeClass(gr,annoDir)

# annotate with overlapping features
gr <- featureAnn(gr, annoDir)

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
gr <- countingBiases(gr,cl,segLocation)
stopCluster(cl)

#Old methylation function - may still be usefull, picks up a lot more methylation
gr <-methylation1(gr,annoDir)

#New methylation functions
gr <- methylation2(gr,annoDir)
gr <- methylationDiff(grannoDir) #Almost no results - very small datasets

#Extra Current annotations
gr <- strainSpec(gr, loci)
gr <- lifeCycle(gr, loci)
gr <- mutantSpec(gr, loci)
gr <- phaseMatch(gr)

#Other functions which were done for arabidopsis

#gr <- histoneAnnotate(gr)
#cl <- makeCluster(24)
#gr <- methAnnotate(gr, cl)
#stopCluster(cl)
#gr <- tissueSpec(gr, loci) #adapted for zygotes, key mutants, strains
#gr <- agoIP(gr, loci)
#cl <- makeCluster(24)
#Pol45(gr, cl)
#stopCluster(cl)
#gr <- annPol(gr)
#Save file
save(gr, loci, file = file.path(segLocation, paste0("gr_fdr", fdr, "_", gitfingerprint, ".RData")))
