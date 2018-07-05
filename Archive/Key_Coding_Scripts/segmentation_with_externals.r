### Mass segementation of all available Chlamy data
#Author: nem34 (Using code from bacms2 previous segmentation)
#Date: 12/11/15


#Import Segment Seq
library(segmentSeq)

#Define data directory
datadir <- "/projects/nick_matthews/segmentation_2018"


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

# as of Adrians mail:

#- I would not use those libraries made by 454 technology. Nobody is using this technology since, I would say, 8 years ago. It is not a bad technology, but it has quite strong bias to certain sequences so that comparison with solexa/illumina is imposible. I would forget about them.
#
#- Please, note that only Betty and I used the strain CC1883. As far as I remember, Andrew Bassett and Attila Molnar used cw15-302 (CC4350 if you want to use the code from Chlamy.org). The info about Daissy's library is correct but need further clarifications (see below).
#
#- Libraries that I run in Eric Miska laboratory are already present in your list! It means that I actually uploaded them into the DCB pipeline (SL2310-to-SL2333).
#
#- I would not use the libraries SL2362-to-SL2373. I made them, and because of the weird treatments that I did to the RNA before generating the libraries, the quality resulted pretty bad. 
#However, in the case that you want to go ahead with them, please note that description of the treatment is not correct!!! -> "Ant Phosp ( ) PNK ( )" Please, have a look at the attached excel file for clarification of the real treatments.
#Regarding biological replicates of these libraries (SL2362-to-SL2373), let say that there are no such a thing in this case. 
#The same sample (either CC1883 or mut47) was treated in different manner with Antarctic phosphatase (Ant Phosp) and/or Polynucleotide Kinase (PNK), but I included no replicates for each treatment. For more info, take a look at the enclosed excel file.
#
#- Indeed, samples SL2126, 2127 and 2128 are biological replicates of the same strain (RC_L2-1 J - cc125). The name means it is the result of a cross (RC means “recombinant”) between the wild type “Japanese” strain (J) and the common wild type CC125). 
#Please, note that ALL the libraries with a “RC” tag in their name (SL2126-to-SL2137) are “RC”, meaning they are all recombinant (crosses) between the Japanese strain and the CC125 (this information is MISSED in the table that you sent me!).
#
#- Even though you now can easily see which library comes from ago3 mutants in the table that I am sending here, I let you know that m33, m49, and m25 correspond to different ago3 mutant alleles.
#

library(readr)
library(dplyr)
meta <- read_csv("/projects/nick_matthews/chlamy_locus_map_github/Summary_of_Data.csv")

## taking out SL2362-73 and SL2126-2137
# they are anntated in the csv file
meta <- meta %>% 
  filter(InCurrentLociRun == "Yes")

#creating alignment objext
aD <- readGeneric(files = meta$File, dir = "",replicates = meta$Replicate, libnames = meta$DataCode,chrs = chrs, chrlens = chrlens,cols=cols, verbose=TRUE, gap = 200,cl=cl)
save(aD, file=file.path(datadir ,"aD_chlamy_segmentation.RData"))
#load("aD_first_chlamy_segmentation_nick.RData")

#Get rid of highly expressed and repetative data
fivenum( aD@alignments$multireads )
# [1]   1   5  35 179 994
aD<-aD[aD@alignments$multireads<100]
#load("aDlt20_first_chlamy_segmentation_nick.RData")

#Process alignment data to find potential segements
sD<-processAD(aD,gap=100, cl=cl) #How big a gap should I use?
save(sD, file=file.path(datadir ,"sD_chlamy_segmentation.RData"))
#load("sDlt20_first_chlamy_segmentation_nick.RData")

#Generate Locus map
hS<-heuristicSeg(sD=sD,aD=aD,getLikes=TRUE,cl=cl)
save(hS, file=file.path(datadir ,"hS_chlamy_segmentation.RData"))

#Generate a genome map
segD<-classifySeg(aD=aD,sD=sD,cD=hS, getLikes=TRUE,cl=cl)
save(segD, file=file.path(datadir ,"segD_chlamy_segmentation.RData"))
#load("segD_first_chlamy_segmentation_nick.RData")
