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
load(file.path(annoDir, annoFile))
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
                     RepClass=ordered(cut(aD@alignments$multireads, breaks=c(0,1,2,5,10,20,50,Inf)))) %>%
  mutate(AbundanceClass = cut(value, breaks=c(0,5,500, 1e3,1e5, 1e6, Inf)))
  return(mydf)
}

# wt <- classSegLike@annotation$Plant.Expt.Type %in% c("WT","Col/Col","dpi10")
ath$matches <- ath$aD@alignments$matches
# select only wt libraries
ath$aDnormal <- ath$aD[, ath$wt]
# correct for double counting of the same reads to multiple positions
ath$aD_wt_dedup <- ath$aDnormal[!duplicated(as.character(ath$aDnormal@alignments$tag)),]
ath$aD_wt_dedup@alignments$multireads <- ath$aD_wt_dedup@alignments$matches

# load chlamy sequences
wt_chlamy <- meta$Controls == "wt" # only WT libs
load(file.path(segLocation, aDfile), envir = chlamy <- new.env())
chlamy$aDnormal <- chlamy$aD[, wt_chlamy]
chlamy$tmp <- chlamy$aDnormal[!duplicated(as.character(chlamy$aDnormal@alignments$tag)),]
abudant_srna <- "TTAGTGACGCGCATGAA" ##abudant culprit!!
chlamy$aD_wt_dedup <- chlamy$tmp[!as.character(chlamy$tmp@alignments$tag) == abudant_srna,]


# how many totat sRNAs? (redundant and non-redundant)
tmp <- chlamy$aD[!duplicated(as.character(chlamy$aD@alignments$tag)),]
sum(tmp@data) / 1e6
nrow(tmp@data) / 1e6
# [1] 9.684103
sum(tmp@data) / 1e6
# [1] 336.6617
# Arabidopsis
tmp2 <- ath$aD[!duplicated(as.character(ath$aD@alignments$tag)),]
sum(tmp2@data) / 1e6
> sum(tmp2@data) / 1e6
[1] 264.3087
> nrow(tmp2@data) / 1e6
[1] 16.29592

# investigating 17bp issue
# chlamy$aDnormal <- chlamy$aD[, "SL2108"]
# identify high 17 species in one lib:
# aDnormal@data <-  aDnormal@data/aDnormal@alignments$multireads
# WT?
# Slot "libnames":
#  [1] "SL2108" "SL2121" "SL2122" "SL2123" "SL2124" "SL2125"
#  [7] "SL2181" "SL2182" "SL2183" "SL2184" "SL2185" "SL2186"
# [13] "SL2187" "SL2188" "SL2189" "SL2301" "SL2302" "SL2303"
# [19] "SL2310" "SL2311" "SL2312" "SL2313" "SL2314" "SL2315"
# [25] "SL2322" "SL2323" "SL2324" "SL2325" "SL2326" "SL2327"
# ok 2108-2125, 2310-23
# high in 17bp: "SL2181-89 2213..  2301.2309
chlamy$aDnormal <- chlamy$aD[, "SL2181"]
tmp <- chlamy$aDnormal[width(chlamy$aDnormal@alignments$tag)==17, ]
myord <- order(tmp@data[,1], decreasing = T)
tmp2 <- tmp[myord,]
myseq <- "TTAGTGACGCGCATGAA" ##abudant culprit!!
myseq <- "TAGTGACGCGCATGAA" # high 17bp
myseq <- "GAATTAAGGCGTACGG" # not in genome
myseq <- "TCGAACCATCTAGTAGCTG" #8
chlamy$aDnormal[as.character(chlamy$aDnormal@alignments$tag)==myseq, ]
chlamy$aD[as.character(chlamy$aD@alignments$tag)==myseq, ]
chlamy$aD_wt_dedup[as.character(chlamy$aD_wt_dedup@alignments$tag)==myseq, ]
# grep "^TCGAACCATCTAGTAGCTG$" SL2108_L1.RUN540.trimmed_reads.fastq | wc -l #8
# finish identify high 17


df_ath    <- aD2df(ath$aD_wt_dedup)
df_chlamy <- aD2df(chlamy$aD_wt_dedup)

round(prop.table(table(width(chlamy$aD_wt_dedup@alignments))[1:15]),3)*100
# non-redundant
   # 15    16    17    18    19    20    21    22    23    24    25    26    27    28
 # 2.66  2.85  3.08  2.75  3.79 35.82 36.53  3.21  2.21  1.79  1.51  1.33  1.19  1.00


dfrep_ath <- df_ath %>%
  group_by(Size, RepClass) %>%
  summarise(non_redundant = sum(value > 0), 
            redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(non_redundant_pct=non_redundant/sum(non_redundant)) %>%
  mutate(Plant = "A. thaliana")

dfNuc_ath <- df_ath %>%
  filter(!is.na(AbundanceClass)) %>%
  group_by(Size, firstNuc) %>%
  summarise(non_redundant = sum(value > 0), 
            redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(non_redundant_pct=non_redundant/sum(non_redundant)) %>%
  mutate(Plant = "A. thaliana")

dfAbund_ath <- df_ath %>%
  filter(!is.na(AbundanceClass)) %>%
  group_by(Size, AbundanceClass) %>%
  summarise(non_redundant = sum(value > 0), 
            redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(non_redundant_pct=non_redundant/sum(non_redundant)) %>%
  mutate(Plant = "A. thaliana")

dfrep_chlamy <- df_chlamy %>%
  group_by(Size, RepClass) %>%
  summarise(redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(Plant = "C. reinhardtii")

dfNuc_chlamy <- df_chlamy %>%
  group_by(Size, firstNuc) %>%
  summarise(redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(Plant = "C. reinhardtii")

dfAbund_chlamy <- df_chlamy %>%
  filter(!is.na(AbundanceClass)) %>%
  group_by(Size, AbundanceClass) %>%
  summarise(non_redundant = sum(value > 0), 
            redundant = sum(value)) %>%
  mutate(redundant_pct=redundant/sum(redundant)) %>%
  mutate(non_redundant_pct=non_redundant/sum(non_redundant)) %>%
  mutate(Plant = "C. reinhardtii")

dfrep <- rbind(dfrep_chlamy, dfrep_ath) %>%
  filter(Size < 30) %>%
  mutate(Plant = as_factor(Plant))

dfNuc <- rbind(dfNuc_chlamy, dfNuc_ath) %>%
  filter(Size < 30) %>%
  mutate(Plant = as_factor(Plant))

dfAbund <- rbind(dfAbund_chlamy, dfAbund_ath) %>%
  filter(Size < 30) %>%
  mutate(Plant = as_factor(Plant))

# dfmelt <- mydf %>%
#   group_by(Size) %>%
#   summarise(redundant = sum(value)/1e6)

ggsize <- ggplot(dfAbund, aes(x=factor(Size),
                        fill=AbundanceClass,
                        y=redundant)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  # xlab("sRNA size") +
  ylab("sRNA read count") +
  theme_bw() +
  # theme(legend.direction="horizontal") +
  theme(legend.justification = c(0, 0), 
        legend.position = c(0, 0),
        axis.title.x=element_blank()) +
  # ggtitle('Size distribution of redundant sRNAs reads') +
  guides(fill = guide_legend(title = "repeat class")) +
  scale_y_continuous(labels = scales::unit_format(unit = "M", scale = 1e-6, digits = 2), 
                     breaks = scales::pretty_breaks(n = 6))
 # ggsave(ggsize,file=file.path(inputLocation,"both_sRNA_redundant_distributionvsmultimatching.pdf"),width=5,height=4)

ggAbund <- ggplot(dfAbund, aes(x=factor(Size),
                        fill=AbundanceClass,
                        y=redundant_pct)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  # xlab("sRNA size") +
  ylab("Abundances of sRNA species") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  # theme(legend.direction="horizontal") +
  theme(legend.position="none",
        axis.title.x=element_blank()) +
  # ggtitle('Size distribution of redundant sRNAs reads') +
  guides(fill = guide_legend(title = "repeat class"))
 # ggsave(ggsize,file=file.path(inputLocation,"both_sRNA_redundant_distributionvsmultimatching.pdf"),width=5,height=4)

ggrep <- ggplot(dfrep, aes(x=factor(Size),
                        fill=RepClass,
                        y=redundant_pct)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  # xlab("sRNA size") +
  ylab("Repetitiveness classes proportions") +
  theme(legend.justification = c(0, 0), 
        legend.position = c(0, 0),
        axis.title.x=element_blank()) +
  theme_bw() +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0)) +
  scale_y_continuous(labels = scales::percent)
# ggsave(gg,file=file.path(inputLocation,"both_sRNA_redundant_distributionvsmultimatching_proportions.pdf"),width=5,height=4)

ggfirstNuc <- ggplot(dfNuc, aes(x=factor(Size),
                        fill=firstNuc,
                        y=redundant_pct)) +
  geom_bar(stat="identity") +
  facet_grid(Plant ~ ., scale="free") +
  scale_fill_manual(values=rev(brewer.pal(7,"RdYlBu"))) +
  scale_fill_brewer(palette = "Paired") +
  xlab("sRNA size") +
  ylab("5' Nucleotide proportions") +
  theme_bw() +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.justification = c(0, 0), legend.position = c(0, 0)) +
  guides(fill = guide_legend(title = "Nucleotide"))
# ggsave(gg,file=file.path(inputLocation,"both_sRNA_redundant_size_Nuc_proportions.pdf"),width=5,height=4)

gg <- grid.arrange(ggsize, ggAbund, ggrep, ggfirstNuc)
ggsave(gg,file=file.path(inputLocation, "sRNA_size_rep_firstNuc_chlamy_vs_arabidopsis.pdf"),width=7,height=7)

table(firstNucnomulti)/length(firstNucnomulti)
#    A    C    G    T 
# 0.38 0.14 0.19 0.28
 round(tapply(rowSums(aDnormalnomulti@data),firstNucnomulti,sum)/sum(aDnormalnomulti@data),2)
#    A    C    G    T 
# 0.37 0.13 0.23 0.28

# do percentage!
df_chlamy %>%
  filter(Size < 30) %>%
  group_by(Size) %>%
  summarise(Count_nonred = sum(value==1), 
            Count = sum(value)) %>%
  mutate(Count_pct = 100 * Count / sum(Count)) %>%
  mutate(Count_nonred_pct = 100 * Count_nonred / sum(Count_nonred))
 # 3    17  1468870     3.84
 # 4    18  1665041     4.35
 # 5    19  2151513     5.62
 # 6    20  3089187     8.07
 # 7    21 10691254    27.9
 # 8    22  3615946     9.44
 # 9    23  1653027     4.32
# 10    24   995817     2.60

# how many redundant 5' nuc sRNAs (%):
df_ath %>%
  filter(Size < 30) %>%
  group_by(RepClass) %>%
  summarise(Count_nonred = sum(value==1), 
            Count = sum(value)) %>%
  mutate(Count_pct = 100 * Count / sum(Count)) %>%
  mutate(Count_nonred_pct = 100 * Count_nonred / sum(Count_nonred))
  # RepClass Count_nonred    Count Count_pct Count_nonred_pct
# 1 (0,1]          231331 16398180    42.8              70.2
# 2 (1,2]           32650  6776050    17.7               9.90
# 3 (2,5]           32766 11946240    31.2               9.94
# 4 (5,10]          11467  1144325     2.99              3.48
# 5 (10,20]          7846  1620860     4.23              2.38
# 6 (20,50]          8652   281034     0.734             2.62
# 7 (50,Inf]         4931   132770     0.347             1.50

# how many redundant 5' nuc sRNAs (%):
df_ath %>%
  filter(Size < 30) %>%
  group_by(firstNuc) %>%
  summarise(Count_nonred = sum(value==1), 
            Count = sum(value)) %>%
  mutate(Count_pct = 100 * Count / sum(Count)) %>%
  mutate(Count_nonred_pct = 100 * Count_nonred / sum(Count_nonred))
  # firstNuc Count_nonred    Count Count_pct Count_nonred_pct
  # <fct>           <int>    <dbl>     <dbl>            <dbl>
# 1 A             2269101 53388019     36.2              38.4
# 2 C              821875 13809555      9.35             13.9
# 3 G             1150109 29312473     19.8              19.4
# 4 T             1675056 51163090     34.6              28.3


df_chlamy %>%
  filter(Size < 30) %>%
  group_by(AbundanceClass) %>%
  summarise(Count_nonred = sum(value==1), 
            Count = sum(value)) %>%
  mutate(Count_pct = 100 * Count / sum(Count)) %>%
  mutate(Count_nonred_pct = 100 * Count_nonred / sum(Count_nonred))
  # AbundanceClass Count_nonred    Count Count_pct Count_nonred_pct
  # <fct>                 <int>    <dbl>     <dbl>            <dbl>
# 1 (0,5]                329643  1425244      3.72              100
# 2 (5,500]                   0  5086403     13.3                 0
# 3 (500,1e+03]               0  1529285      3.99                0
# 4 (1e+03,1e+05]             0 17925920     46.8                 0
# 5 (1e+05,1e+06]             0 10434365     27.2                 0
# 6 (1e+06,Inf]               0  1898242      4.96                0
# 7 NA                        0        0      0                   0

# finding miRNAs in datasets
miRNA_vs_aD <- findOverlaps(miRNA, chlamy$aD_wt_dedup@alignments)
tmp <- colSums(chlamy$aD_wt_dedup[subjectHits(miRNA_vs_aD),]@data)
tmp2 <- colSums(chlamy$aD_wt_dedup@data)

 round(tmp/tmp2 * 100, 2)

# SL2108 SL2121 SL2122 SL2123 SL2124 SL2125 SL2181 SL2182 SL2183 SL2184 SL2185 SL2186 SL2187 SL2188 SL2189 SL2301 SL2302 SL2303 SL2310 SL2311 SL2312 SL2313 SL2314 SL2315
  # 3.95   4.10   3.99   1.02   1.04   1.17   0.88   0.86   1.00   0.79   0.94   0.87   0.68   0.52   0.59   0.75   0.59   0.64   4.54   6.75   8.04   4.23   3.49   4.47
# SL2322 SL2323 SL2324 SL2325 SL2326 SL2327
  # 7.99   7.54   6.82   7.21   8.09   6.30
# miRNA percentage ranges from 0.5 to 8%

range(width(gr))
# [1]    16 15201
cov <- sum(width(gr))
cov/1e6
# [1] 4.571925
genome_length <- sum(width(segD@coordinates))
(cov/genome_length)*100
# [1] 4.128374
genome_length/1e6


