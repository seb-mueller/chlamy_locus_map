###################################################################
# date: 2018-10-18
# chlamy small RNA loci R code
# This code is based on a finished loci definition as obtained from:
# chlamy_annotation_pipeline.R and Segmentation_Analysis.R


library(segmentSeq)
#Specify location of the segmentation
#setwd("C:/Users/Nick/Documents/Uni Work/Third Year/Project/chlamy_locus_map")
#Specify location of segmentation
baseDir <- "/projects/nick_matthews"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
#set working directory to github repository on cluster
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")
#segLocation <- "C:/Users/Nick/Documents/Uni Work/Third Year/Project/segmentation_with_externals.r_2015-11-12_17/"
annoDir <- "/projects/nick_matthews/resources"
#annoDir = "C:/Users/Nick/Documents/Uni Work/Third Year/Project/segmentMap_II/Old_Annotation_Files"
inputdata <- "segD_chlamy_segmentation_multi200_gap100.RData" #13739
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
# prefix  <-  "multi200_gap100"
saveLocation <- file.path(segLocation, prefix)

# code for locus summary plots
# the locus map (not having selected loci yet, but after calculating loci likelihoods) is 'segD'
load(file.path(segLocation, inputdata))
# imports segD object, containing all loci plus count data, posteriours etc for each loci
segD@coordinates
# GRanges object with 13739 ranges and 0 metadata columns:
#               seqnames         ranges strand
#                  <Rle>      <IRanges>  <Rle>
#       [1] chromosome_1 [    1,  3542]      *
#   [13738]  scaffold_50 [11355, 11743]      *
#   [13739]  scaffold_52 [    1,   663]      *
#   -------
#   seqinfo: 54 sequences from an unspecified genome; no seqlengths
load(file.path(saveLocation, "gr_fdr0.05_59a4cbc.RData"))
# imports gr, loci

# check parameter used to generate gr
attr(gr, "parameter")
# [1] "_multi200_gap100"
attr(gr, "git")
# [1] "59a4cbc"
