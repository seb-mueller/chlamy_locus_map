###################################################################
# date: 2018-10-18
# chlamy small RNA loci R code
# This code is based on a finished loci definition as obtained from:
# chlamy_annotation_pipeline.R and Segmentation_Analysis.R
# used to investigate loci for interesting association etc and descriptive features for writing paper

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
library(purrr)

baseDir <- "/projects/nick_matthews"
# baseDir <- "/home/sm934/workspace/chlamy"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
annoFile    <- "chlamy_all_annotations.Rdata"
#set working directory to github repository on cluster
gitdir  <- file.path(baseDir, "chlamy_locus_map_github")
# list.files(segLocation, pattern = "segD.*")
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
aDfile    <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"

gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
# gitfingerprint <- "9dd50d8"
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
# [1] "LociRun2018_multi200_gap100"
saveLocation <- file.path(segLocation, paste(prefix, gitfingerprint, sep = "_"))
dir.create(saveLocation)


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
table2 <- function(...) table(..., useNA = "always")  # counts all three

table2d <- function(...) ftable(addmargins(table(..., useNA = "always")))
#------------------------------------ examining individual features
table2(gr$sizeclass)
#
#      (0,100]     (100,400]    (400,1500] (1500,3e+03] (3e+03,Inf]
#          272          2113          3109         552         118
table2(gr$predominant_sRNA_sizeClass)
#   equal_20bp   equal_21bp  larger_21bp smaller_20bp
#          415         3302         1080         1367
table2(gr$ratio_strand_class)
# strong_bias    med_bias     no_bias        <NA> 
#        2536        2182         829         617 
table2(gr$repetitivenessClass)
# 
#  low  med high <NA> 
#  518 1777 3477  392 
table2(gr$phaseClass)
#   none median   high 
#   6020    123     21 

#------------------------------------ interessing associations

table(gr$predominant_sRNA_sizeClass, gr$repetitivenessClass)
table2d(gr$predominant_sRNA_sizeClass, gr$repetitivenessClass)
#                low  med high  Sum
#                                  
# equal_20bp      12   98  298  408
# equal_21bp     139  780 2368 3287
# larger_21bp    261  441  340 1042
# smaller_20bp   106  458  471 1035
# Sum            518 1777 3477 5772
# -> 20/21 are very Repetitive, but smalle/bigger not
chisq.test(table(gr$predominant_sRNA_sizeClass, gr$predominant_5prime_letter))
ftable(addmargins(table(gr$predominant_sRNA_sizeClass, gr$sizeclass)))
#               (0,100] (100,400] (400,1.5e+03] (1.5e+03,3e+03] (3e+03,Inf]  Sum
#
# equal_20bp         28       187           186              13           1  415
# equal_21bp        144      1106          1606             349          97 3302
# larger_21bp        30       358           565             117          10 1080
# smaller_20bp       70       462           752              73          10 1367
# Sum               272      2113          3109             552         118 6164
ftable(addmargins(table(gr$predominant_sRNA_sizeClass, gr$predominant_5prime_letter)))
#                  A   AC  ACG   AG   AT    C   CG  CGT   CT    G   GT    T  Sum
# equal_20bp      10    2    0    2    3   40    1    0   10  132   18  191  409
# equal_21bp     204   11    0   23  146  125   47   15  114   53  113 2347 3198
# larger_21bp    100   15    5   36   13   99  168    1   19  211   32  216  915
# smaller_20bp   113   17    3   31    3  124  201    1   21  172   24  220  930
# Sum            427   45    8   92  165  388  417   17  164  568  187 2974 5452
# -> most 21bp start with T (no G!), 20bp start T or G
