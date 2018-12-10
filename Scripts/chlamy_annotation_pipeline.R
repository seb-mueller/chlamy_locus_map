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
library(MKmisc) #binomCI, binomial confidence interval
library(readr)
library(dplyr)
# library(simpleboot)

###### setting variables

fdr <- 0.05
baseDir <- "/projects/nick_matthews"
# baseDir <- "/home/sm934/workspace/chlamy"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
annoFile    <- "chlamy_all_annotations.Rdata"
#set working directory to github repository on cluster
gitdir  <- file.path(baseDir, "chlamy_locus_map_github")
meta    <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
  filter(InCurrentLociRun == "Yes")
metawt <- meta$Controls %in% c("wt") # only WT libs

source(file.path(gitdir, "Scripts/chlamy_source_code.R"))

# list.files(segLocation, pattern = "segD.*")
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
aDfile    <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"
# inputdata <- "segD_chlamy_segmentation_smallset_200.RData" # 2loci
# inputdata <- "segD_chlamy_segmentationmulti200_wt_adrian100.RData" #8977
# inputdata <- "segD_chlamy_segmentationmulti200_wt_adrian200.RData" #8158
# using instead of arbitray versions
# e.g. "1f6085a"
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
gitfingerprint <- "8ab6f64"
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
# [1] "LociRun2018_multi200_gap100"
saveLocation <- file.path(segLocation, paste(prefix, gitfingerprint, sep = "_"))
dir.create(saveLocation)


cl <- makeCluster(24)

#####Load in and process segmentation, selecting significant loci#####
# load(file.path(segLocation,"segD_first_chlamy_segmentation_nick.RData"))
# loads segD object (lociData class)
load(file.path(segLocation, inputdata))

# Select loci based on some fdr
# perReplicate: If TRUE, selection of loci is done on a replicate by replicate basis. If FALSE, selection will be done on the likelihood that the locus represents a true locus in at least one replicate group.

# which setting are we using? see issue #9
loci <- selectLoci(cD = segD, FDR = fdr, perReplicate = TRUE) # 6164
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


#####These next few functions compute and compile annotation files, this shouldn't need runnding every time#####
#Compute introns - this takes a while, don't run unless necessary
#intronCalculate()
#Process transposon file

annotenv <- new.env()
myfile <- file.path(annoDir, "transposon_annotations.Rdata")
if(!file.exists(myfile)) {
  transposonProcess(annoDir, gitdir)
} else {
  message(myfile, " exists, loading..")
  load(myfile, envir = annotenv, verbose = FALSE)
}

myfile <- file.path(annoDir, annoFile)
if(!file.exists(myfile)) {
  #Compute and compiles annotations
  compileAnnotations(annoDir, annoFile = annoFile)
} else {
  message(myfile, " exists, loading..")
  load(myfile, envir = annotenv, verbose = FALSE)
}
# attach annotenv so objects are accesible
attach(annotenv)
ls(annotenv)
#  [1] "anno"           "CDS"            "exons"          "fiveprimeUTR"  
#  [5] "genes"          "introns"        "irs"            "miRNA"         
#  [9] "mRNA"           "MSAT"           "promoter"       "repetativeSeq" 
# [13] "rRNA"           "TE_Class_DNA"   "TE_Class_RET"   "TE_Copia"      
# [17] "TE_DIRS"        "TE_DualenRandI" "TE_EnSpm"       "TE_Gulliver"   
# [21] "TE_Gypsy"       "TE_Harbinger"   "TE_hAT"         "TE_L1"         
# [25] "TE_Mariner"     "TE_Novosib"     "TE_Order_DIRS"  "TE_Order_LINE" 
# [29] "TE_Order_LTR"   "TE_Order_SINE"  "TE_Order_TIR"   "TE_P"          
# [33] "TE_REM1"        "TE_RTE"         "TE_TCR1"        "TE_TOC1"       
# [37] "TE_TOC2"        "TE_Undefined"   "threeprimeUTR"  "transposons"   
# [41] "trs"           

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
gr <- featureAnn(locAnn = gr, annotations = annotenv)
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
# TODO: PhaseTank
# gr <- phaseMatch2(gr,annoDir=annoDir,outputName="Pred_tab_2018.11.27_18.08")

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
gr <- countingBiases(locAnn = gr, cl = cl,
                     segLocation = segLocation,
                     wt = metawt,
                     aDfile = aDfile)

# classing continues features using thresholds obtained by inspections
  #For 21 vs 20 ratio
  locAnn$ratio21vs20Class <- classCI(locAnn$counts21, locAnn$counts20, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="21vs20_0.3_0.7.pdf")
  #For 20 vs 21 ratio
  locAnn$ratio20vs21Class <- classCI(locAnn$counts20, locAnn$counts21, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="20vs21_0.3_0.7.pdf")
  #For small RNAs vs Normal RNAs
  locAnn$ratioSmallvsNormalClass <- classCI(locAnn$countsSmall, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="SmallvsNormal_0.3_0.7.pdf")
  #For big RNAs vs Normal RNAs
  locAnn$ratioBigvsNormalClass <- classCI(locAnn$countsBig, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="BigvsNormal_0.3_0.7.pdf")
  # confidence intervals on strands, as before
  locAnn$ratio_strand_class <- classCI(locAnn$countsplus, locAnn$countsminus,
                                       probs = c(0.1,0.35,0.65,0.9), comma = 1,plotname="plusvsminus2.pdf") #changed probs to 0.1 and 0.9 rather than 0.2 and 0.8
  levels(locAnn$ratio_strand_class) <- c("strong bias", "med bias", "no bias", "med bias", "strong bias")

  locAnn$repetitivenessClass <- classCI(rowSums(countsnormalwt) - rowSums(countsnormalwtnorm), rowSums(countsnormalwtnorm), probs = c(med = 0.3, high = 0.7), comma = 1, main = "Repetitiveness")

stopCluster(cl)


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

# export as gff3 file for viewing in browser
export.gff3(gr, con = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".gff")))
#Write csv for phasing

write.csv(as.data.frame(gr), file = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".csv")))
#Save file
save(gr, loci, baseDir, prefix, saveLocation, file = file.path(saveLocation, paste0("gr_fdr", fdr,  ".RData")))
load(file = file.path(saveLocation, paste0("gr_fdr", fdr,  ".RData")))
