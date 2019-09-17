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
library(ggplot2)

baseDir <- "/projects/nick_matthews"
# baseDir <- "/home/sm934/workspace/chlamy"
# Specify location of segmentation
segLocation <- file.path(baseDir, "segmentation_2018")
annoDir     <- file.path(baseDir, "resources")
annoFile    <- "chlamy_all_annotations.Rdata"
lociRun <- "LociRun2018_multi200_gap100_90c7213"
lociLocation <- file.path(baseDir, "segmentation_2018", lociRun)
#set working directory to github repository on cluster
gitdir  <- file.path(baseDir, "chlamy_locus_map_github")
# list.files(segLocation, pattern = "segD.*")
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
aDfile    <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"

gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
# gitfingerprint <- "90c7213"
prefix <- str_replace(inputdata, "segD_chlamy_segmentation_(.*).RData", "\\1")
# [1] "LociRun2018_multi200_gap100"
saveLocation <- file.path(segLocation, paste(prefix, gitfingerprint, sep = "_"))
inputLocation <- file.path(baseDir, "segmentation_2018", lociRun)
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
load(file.path(saveLocation, "gr_fdr0.05.RData"))
wt_chlamy <- meta$Controls %in% c("wt") # only WT libs
# imports gr, loci, meta

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

# comparing loci definition with arabidopsis
load(file.path(baseDir, "resources/arabidopsis/aD.RData"), envir = ath <- new.env())
#redundant: 
load(file.path(baseDir, "resources/arabidopsis/gr9_clustered_update.rdata"), envir = ath )
load(file.path(baseDir, "resources/arabidopsis/lociIV_fdr01.rdata"), envir = ath )
load(file.path(baseDir, "resources/arabidopsis/workspace.rdata"), envir = ath )

# comparing sRNA sizes and 5' (redundant/non redundant)

aD2df <- function(aD) {
  # sum(aD[ath$matches==1, ])/sum(aD)
  # aD@data <-  aD@data/matches
  # tmp <- table(width(aD@alignments))
  mydf <- data.frame(value=rowSums(aD@data),
                     firstNuc=substr(aD@alignments$tag,1,1),
                     Size=width(aD@alignments), 
                     RepClass=ordered(cut(aD@alignments$multireads, breaks=c(0,1,2,5,10,20,50,Inf))))
  return(mydf)
}
# wt <- classSegLike@annotation$Plant.Expt.Type %in% c("WT","Col/Col","dpi10")
ath$matches <- ath$aD@alignments$matches
# tmp <- aDnormal[as.character(ath$aD@alignments$tag)=="GGGTTTAGGGTTTAGGG", ath$wt]
# select only wt libraries
ath$aDnormal <- ath$aD[, ath$wt]
# correct for double counting of the same reads to multiple positions
ath$aD_wt_dedup <- ath$aDnormal[!duplicated(as.character(ath$aDnormal@alignments$tag)),]
ath$aD_wt_dedup@alignments$multireads <- ath$aD_wt_dedup@alignments$matches

# load chlamy sequences
load(file.path(segLocation, aDfile), envir = chlamy <- new.env())
chlamy$aDnormal <- chlamy$aD[, wt_chlamy]
chlamy$aD_wt_dedup <- chlamy$aDnormal[!duplicated(as.character(chlamy$aDnormal@alignments$tag)),]
# aDnormal@data <-  aDnormal@data/aDnormal@alignments$multireads
# WT?
# Slot "libnames":
#  [1] "SL2108" "SL2121" "SL2122" "SL2123" "SL2124" "SL2125"
#  [7] "SL2181" "SL2182" "SL2183" "SL2184" "SL2185" "SL2186"
# [13] "SL2187" "SL2188" "SL2189" "SL2301" "SL2302" "SL2303"
# [19] "SL2310" "SL2311" "SL2312" "SL2313" "SL2314" "SL2315"
# [25] "SL2322" "SL2323" "SL2324" "SL2325" "SL2326" "SL2327"

df_ath <- aD2df(ath$aD_wt_dedup)
df_chlamy <- aD2df(chlamy$aD_wt_dedup)

# do percentage!
dfsize_ath <- df_ath %>%
  group_by(Size) %>%
  summarise(Count = sum(value)/1e6)

dfrep_ath <- df_ath %>%
  group_by(Size, RepClass) %>%
  summarise(Count = sum(value)) %>%
  mutate(Proportion=Count/sum(Count)) %>%
  mutate(Plant = "A. thaliana")

dfNuc_ath <- df_ath %>%
  group_by(Size, firstNuc) %>%
  summarise(Count = sum(value)) %>%
  mutate(Proportion=Count/sum(Count)) %>%
  mutate(Plant = "A. thaliana")

dfrep_chlamy <- df_chlamy %>%
  group_by(Size, RepClass) %>%
  summarise(Count = sum(value)) %>%
  mutate(Proportion=Count/sum(Count)) %>%
  mutate(Plant = "C. reinhardtii")

dfNuc_chlamy <- df_chlamy %>%
  group_by(Size, firstNuc) %>%
  summarise(Count = sum(value)) %>%
  mutate(Proportion=Count/sum(Count)) %>%
  mutate(Plant = "C. reinhardtii")

dfrep <- rbind(dfrep_chlamy, dfrep_ath) %>%
  filter(Size < 30)

dfNuc <- rbind(dfNuc_chlamy, dfNuc_ath) %>%
  filter(Size < 30)

# dfmelt <- mydf %>%
#   group_by(Size) %>%
#   summarise(Count = sum(value)/1e6)

gg <- ggplot(dfrep, aes(x=factor(Size),
                        fill=RepClass,
                        y=Count)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  xlab("sRNA size") +
  ylab("sRNA read count [Millions]") +
  theme_bw() +
  # ggtitle('Size distribution of redundant sRNAs reads') +
  guides(fill = guide_legend(title = "repeat class")) +
  scale_y_continuous(labels = scales::unit_format(unit = "", scale = 1e-6, digits = 2), 
                     breaks = scales::pretty_breaks(n = 8))
ggsave(gg,file=file.path(inputLocation,"both_sRNA_redundant_distributionvsmultimatching.pdf"),width=5,height=4)

gg <- ggplot(dfrep, aes(x=factor(Size),
                        fill=RepClass,
                        y=Proportion)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  xlab("sRNA size") +
  ylab("Proportion") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = "repeat class"))
ggsave(gg,file=file.path(inputLocation,"both_sRNA_redundant_distributionvsmultimatching_proportions.pdf"),width=5,height=4)

gg <- ggplot(dfNuc, aes(x=factor(Size),
                        fill=firstNuc,
                        y=Proportion)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  xlab("sRNA size") +
  ylab("Proportion") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  guides(fill = guide_legend(title = "repeat class"))
ggsave(gg,file=file.path(inputLocation,"both_sRNA_redundant_size_Nuc_proportions.pdf"),width=5,height=4)

ggsave(gg,file=file.path(inputLocation,"Clustercoverage.eps"),width=30,height=5)


table(firstNucnomulti)/length(firstNucnomulti)
#    A    C    G    T 
# 0.38 0.14 0.19 0.28
 round(tapply(rowSums(aDnormalnomulti@data),firstNucnomulti,sum)/sum(aDnormalnomulti@data),2)
#    A    C    G    T 
# 0.37 0.13 0.23 0.28
