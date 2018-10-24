### Mass segementation of all available Chlamy data
#Author: nem34, sm934 (Using code from bacms2 previous segmentation)
#Date: 06/07/2018
# how to run:
# in git Script directory (on node7):
# /scripts/csmit -m 100G "/applications/R/R-3.5.0/bin/Rscript chlamy_segmentation_gap100.R"

library(segmentSeq)
library(readr)
library(dplyr)

#Define crucial parameters
# /data/pipeline/prod/
baseDir      <- "/projects/nick_matthews"
segLocation  <- file.path(baseDir, "segmentation_2018")
gitdir       <- file.path(baseDir, "chlamy_locus_map_github")
mymultireads <- 200
mygap        <- 100
mycomment    <- paste0("LociRun2018_multi", mymultireads, "_gap", mygap) # goes in filenames

# selected libraries as discussed with Nick/Seb/Adrian (july2018)
meta <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
    filter(LociRun2018 == "Yes")

# using instead of arbitray versions
# e.g. "1f6085a"
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)

#Make cluster
cl<-makeCluster(46)

#define Chromosome lengths
chrlens <- c(8033585,9223677,9219486,4091191,3500558,9023763,6421821,5033832,7956127,6576019,
		3826814,9730733,5206065,4157777,1922860,7783580,7188315,271631,219038,200793,
		189560,163774,127913,127161,102191,80213,55320,55278,52813,52376,
		48183,42264,39192,33576,32450,25399,24537,24437,22408,22082,
		21325,21000,20974,17736,16939,16627,14746,14165,13462,
		12727,11225,6241,2479,2277)

#Define chromosomes and scaffolds
chrs <- c("chromosome_1","chromosome_2","chromosome_3","chromosome_4","chromosome_5","chromosome_6","chromosome_7","chromosome_8","chromosome_9","chromosome_10",
		"chromosome_11","chromosome_12","chromosome_13","chromosome_14","chromosome_15","chromosome_16","chromosome_17","scaffold_18","scaffold_19","scaffold_20",
		"scaffold_21","scaffold_22","scaffold_23","scaffold_24","scaffold_25","scaffold_26","scaffold_27","scaffold_28","scaffold_29","scaffold_30",
		"scaffold_31","scaffold_32","scaffold_33","scaffold_34","scaffold_35","scaffold_36","scaffold_37","scaffold_38","scaffold_39","scaffold_40",
		"scaffold_41","scaffold_42","scaffold_43","scaffold_44","scaffold_45","scaffold_46","scaffold_47","scaffold_48","scaffold_49","scaffold_50",
		"scaffold_51","scaffold_52","scaffold_53","scaffold_54")

#Define Columns in Input Data
cols = c(chr = 1, tag = 10, start = 4, end = 5, count = 6, strand = 7)

#creating alignment objext
# file.exists(file.path(baseDir, "sequencing_data", meta$File))

 # meta$File[!file.exists(file.path(baseDir, "sequencing_data", meta$File))]

aD <- readGeneric(files = meta$File,
                  dir = file.path(baseDir, "sequencing_data"),
                  replicates = meta$Replicate,
                  libnames = meta$DataCode,chrs = chrs,
                  chrlens = chrlens,
                  cols=cols,
                  verbose=TRUE,
                  gap = mygap,
                  cl=cl)

# saving version information within the R object
attr(aD, "parameter")   <- mycomment
attr(aD, "git")         <- gitfingerprint
attr(aD, "sessionInfo") <- sessionInfo()

#Get rid of highly expressed and repetative data
# fivenum( aD@alignments$multireads )
# [1]   1   5  35 179 994
aD<-aD[aD@alignments$multireads < mymultireads]
save(aD, meta, file=file.path(segLocation ,paste0("aD_chlamy_segmentation_",mycomment,".RData")))
#load("aDlt20_first_chlamy_segmentation_nick.RData")

#Process alignment data to find potential segements
sD<-processAD(aD,gap=mygap, cl=cl) #How big a gap should I use?
 save(sD, meta, file=file.path(segLocation ,paste0("sD_chlamy_segmentation_",mycomment,".RData")))
#load(file=file.path(segLocation ,paste0("sD_chlamy_segmentation_",mycomment,".RData")))

#Generate Locus map
hS<-heuristicSeg(sD=sD,aD=aD,getLikes=TRUE,cl=cl)
save(hS, meta, file=file.path(segLocation, paste0("hS_chlamy_segmentation_", mycomment, ".RData")))

#Generate a genome map
segD<-classifySeg(aD=aD,sD=sD,cD=hS, getLikes=TRUE,cl=cl)
attr(segD, "parameter")   <- mycomment
attr(segD, "git")         <- gitfingerprint
attr(segD, "sessionInfo") <- sessionInfo()
save(segD, meta, file=file.path(segLocation ,paste0("segD_chlamy_segmentation_", mycomment, ".RData")))
