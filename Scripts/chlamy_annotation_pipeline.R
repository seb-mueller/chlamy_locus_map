#Script for running annotation functions as well as defining sRNA loci for chlamy.
# The appropriate thresholds (e.g. FDR) are determined in:
# Segmentation_Analysis.R
#Adapted by Nick Matthews from code written by Tom Hardcastle and Sebastian Muller
#Date: 09/02/16
##### Load necesary libraries
library(stringr)
library(xtable)
library(rtracklayer)
library(segmentSeq)
# library(pROC)
library(MASS)
library(RColorBrewer)
# library(mclust)
library(baySeq)
library(MKmisc)
library(readr)
library(dplyr)
# library(simpleboot)

###### setting variables

fdr <- 0.05
baseDir <- "/projects/nick_matthews"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
annoFile <- "chlamy_all_annotations.Rdata"
#set working directory to github repository on cluster
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")
meta <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
    filter(InCurrentLociRun == "Yes")
#TODO we need to check definition of WTs
metawt <- meta$Controls %in% c("wt")

source(file.path(gitdir, "Scripts/chlamy_source_code.R"))

# list.files(segLocation, pattern = "segD.*")
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #13739
aDfile <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"
# inputdata <- "segD_chlamy_segmentation_smallset_200.RData" # 2loci
# inputdata <- "segD_chlamy_segmentationmulti200_wt_adrian100.RData" #8977
# inputdata <- "segD_chlamy_segmentationmulti200_wt_adrian200.RData" #8158
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
# [1] "multi200_gap100"
saveLocation <- file.path(segLocation, prefix)
dir.create(saveLocation)

# using instead of arbitray versions
# e.g. "1f6085a"
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)

cl <- makeCluster(24)

#####Load in and process segmentation, selecting significant loci#####
# load(file.path(segLocation,"segD_first_chlamy_segmentation_nick.RData"))
# loads segD object (lociData class)
load(file.path(segLocation, inputdata))

# Select loci based on some fdr
# perReplicate: If TRUE, selection of loci is done on a replicate by replicate basis. If FALSE, selection will be done on the likelihood that the locus represents a true locus in at least one replicate group.

# which setting are we using? see issue #9
loci <- selectLoci(cD = segD, FDR = fdr, perReplicate = TRUE) # 4915
attr(loci, "parameter") <- prefix
attr(loci, "git") <- gitfingerprint
# loci <- selectLoci(cD = segD, FDR = fdr, perReplicate = TRUE) # 4915 (plus warning)

## creating/exporting coordinates object
gr <- loci@coordinates

# safe attributes in object
attr(gr, "parameter") <- prefix
attr(gr, "git") <- gitfingerprint
# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr) <- sprintf("CRSL%05.f0", 1:length(gr))

#load("gr_just_cb.RData") #load gr you want to add to it...

# export as gff3 file for viewing in browser
export.gff3(gr, con = file.path(segLocation, paste0("loci_fdr", fdr, prefix, ".gff")))
#import.gff3("loci_fdr01.gff3")
#Write csv for phasing
write.csv(as.data.frame(gr),file=file.path(baseDir, "phasing/loci_for_phasing.csv"))

#####These next few functions compute and compile annotation files, this shouldn't need runnding every time#####
#Compute introns - this takes a while, don't run unless necessary
#intronCalculate()
#Process transposon file

transposonProcess(annoDir,gitdir)
#Compute and compiles annotations
compileAnnotations(annoDir, annoFile = annoFile)
#####the next set of functions take the annotated locus object 'gr' and add some more annotation data#####
#The functions are found 'chlamy_source_code.R'

# annotate by size class
gr <- sizeClass(gr, intervals = c(0,100,400,1500,3000,Inf))
table(gr$sizeclass)
#
#      (0,100]     (100,400]    (400,1500] (1500,3e+03] (3e+03,Inf]
#          272          2113          3109         552         118

# annotate with overlapping features
gr <- expressionClass(locAnn = gr, loci = loci, wt = metawt)
gr <- featureAnn(locAnn = gr, annoDir = annoDir, annoFile = annoFile)

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
gr <- countingBiases(locAnn = gr, cl = cl,
                     segLocation = segLocation,
                     wt = metawt,
                     aDfile = aDfile)
stopCluster(cl)

#Old methylation function - may still be usefull, picks up a lot more methylation
# gr <-methylation1(gr,annoDir)

#New methylation functions
# gr <- methylation2(gr,annoDir)
gr <- methylation(gr, annoDir)
#gr <- methylationDiff(gr,annoDir) #Almost no results - very small datasets
#gr <- methylationDiff(grannoDir) #Almost no results - very small datasets

#Extra Current annotations
gr <- strainSpec(gr, loci, meta = meta, gitdir = gitdir)
gr <- lifeCycle(gr, loci, meta = meta, gitdir = gitdir)
gr <- mutantSpec(gr, loci, meta = meta, gitdir = gitdir)
# gr <- phaseMatch(gr)

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
save(gr, loci, baseDir, prefix, saveLocation, file = file.path(saveLocation, paste0("gr_fdr", fdr, "_", gitfingerprint,  ".RData")))
