#Code adapted by Nick Matthew from code written by Tom Hardcstle
#Analyses quality of segmentation producing a variety of graphs
#Date: 06/07/18

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

# code for locus summary plots mostly to determine fdr cut-off
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

# check parameter used to generate gr
attr(gr, "parameter")
# [1] "_multi200_gap100"
attr(gr, "git")
# [1] "59a4cbc"

# nct<-segD
                                        
master <- read.csv(file.path(gitdir, "Summary_of_Data.csv"))
#To select apparent 'WT' - i.e. non-mutants
#wt<-subset(master,master[,"Ecotype"] =="wt")
#wtReps<-unique(wt[,"Replicate"])

#To select very specific control WTs
wt<-subset(master,master[,"Controls"] =="wt")
wtReps<-unique(wt[,"Replicate"]) #Currently 10 06/07/18
#  [1] 31 32 39 40 41 55 58 59 62 63 # 15/8/18

nfdrnum <- NULL

# check how many loci get selected at various FDR levels
for(ii in (1:20 / 4)) {
  fdr <- 10^-ii
  loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE) #
  nfdrnum <- cbind(nfdrnum, c(fdr, nrow(loci)))
}

png("fdrNumbers.png")
plot(nfdrnum[1,], nfdrnum[2,], log = "xy", xlab = "FDR", ylab = "# of loci", col = rep(c("red", "blue"), each =20))
dev.off()

# check how many loci appear in N replicate groups (at various fdr levels)
pdf("fdr_hists.pdf", width = 10, paper = "a4r")
par(mfcol = c(3,5))
for(ii in 1:5) {
  fdr <- 10^-ii
  loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE)
  hist(rowSums(exp(loci@locLikelihoods[,wtReps])), breaks = 0:length(wtReps), main = paste("FDR =", fdr), xlab = "Expectation")
  hist(rowSums(exp(loci@locLikelihoods)), breaks = 0:nlevels(segD@replicates), main = paste("FDR =", fdr), xlab = "Expectation")
  # also have a look at locus lenght distribution while we're at it.
  plot(density(log10(width(loci@coordinates))), main = "", xlab = "log width")
}
dev.off()

                                        # pick an FDR and move on

#########nloci <- loci[,-1]; nloci@locLikelihoods <- loci@locLikelihoods[,-1]

fdr <- 0.05
loci <- selectLoci(segD, FDR = fdr, perReplicate = TRUE)

# get an idea of the sequencing depth added by each replicate groups
sumLibsizes <- sapply(levels(loci@replicates), function(rep) sum(libsizes(loci)[loci@replicates == rep]))
meanLibscale <- sapply(levels(loci@replicates), function(rep) mean(libsizes(loci)[loci@replicates == rep]))

# order by increasing sequencing depth added
ordLoc <- order(sumLibsizes, decreasing = FALSE)

# number of additional loci added by each library (with increasing sequencing depth)
cumloc <- sapply(1:length(ordLoc), function(ii) {
  message(ii)
  selLoc <- segD[,segD@replicates %in% ordLoc[1:ii]]
  selLoc@locLikelihoods <- as.matrix(selLoc@locLikelihoods[,ordLoc[1:ii],drop = FALSE])
  z <- try(selectLoci(selLoc, FDR = 0.05 , perReplicate = TRUE))
    if(class(z) == "try-error") return(0) else return(nrow(z))
})

collibs <- rep("black", nlevels(loci@replicates))
collibs[wtReps] <- "red"

# plot number of loci discovered as we add deeper libraries. Ideally, we want to see the number of additional loci tailing off, indicating we've achieved enough sequencing depth/variety to get most loci
png("CumSeqVolume.png", height = 800, width = 600)
plot(x = cumsum(sumLibsizes[ordLoc]), y = cumloc, xlab = "Cumulative sequencing volume", ylab = "Total loci discovered", col = collibs[ordLoc], pch = 19)
dev.off()

# individual numbers of loci per replicate group. Should correlate roughly with library size for lower library sizes, hopefully become more or less constant for higher library sizes as returns diminish.
png("LibScalingFactor.png", height = 800, width = 600)
plot(y = summariseLoci(loci, perReplicate = TRUE), x = meanLibscale, log = "x", ylab = "# of loci", xlab = "Library scaling factor")
dev.off()
