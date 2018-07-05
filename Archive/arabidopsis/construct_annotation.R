setwd("/home/tjh48/Code/segmentMap_III")

# I think most of these libraries aren't necessary anymore. Feel free to try without them and load when you get errors.

library(clv)
library(LabelCompare)
library(xtable)
library(rtracklayer)
library(reshape)
library(segmentSeq)
library(org.At.tair.db)
library(pROC)
library(MASS)
library(RColorBrewer)
library(Kendall)
library(RANN)
library(mclust)
source("/home/sm934/code/R/seb_functions.r")
source("construct_annotation_functions.R")

# Load up segmentation
load("nct_classSegLike.RData")

# Select loci based on some fdr
loci9 <- selectLoci(classSegLike, FDR = 0.05, perReplicate = TRUE) #

## creating/exporting coordinates object
gr9 <- loci9@coordinates 

png("example_binning.png", height = 800, width = 600)
par(mfrow = c(2,1))
plot(density(gr9$repetitiveness, na.rm = TRUE, from = 0, to = 1), xlab = "repetitiveness", main = ""); abline(v = c(0.3, 0.6), lty = 2, col = "red")
plot(density(log10(width(gr9))), xlab = "log10 locus width", main = ""); abline(v = log10(c(50, 2000)), lty = 2, col = "red")
dev.off()

# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr9) <- sprintf("CRSL%05.f0",1:length(gr9))

# export as gff3 file for viewing in browser
export.gff3(gr9,con="loci9_fdr01.gff3")

# the next set of functions are found 'construct_annotation_functions.R'. They take the annotated locus object 'gr9' and add some more annotation data. Each of these functions will need looking at and likely tweaking for chlamy.

# annotate by size class
gr9 <- sizeClasses(gr9)

# annotate with overlapping features
gr9 <- featureAnn(gr9, loci9)

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 24s, or a higher number of reads starting with As than usual
cl <- makeCluster(24)
gr9 <- countingBiases(gr9,cl)
stopCluster(cl)

save(gr9, file = "gr9_sc_fa_cb.RData")



# this stuff may not be relevant to the chlamy data, so I haven't annotated it. Let me know when you've got through the stuff above, and we can look at where to go next.


gr9 <- histoneAnnotate(gr9)

gr9 <- phaseMatch(gr9)

cl <- makeCluster(24)
gr9 <- methAnnotate(gr9, cl)
stopCluster(cl)

gr9 <- tissueSpec(gr9, loci9)
gr9 <- agoIP(gr9, loci9)

cl <- makeCluster(24)
Pol45(gr9, cl)
stopCluster(cl)

gr9 <- annPol(gr9)

gr9 <- easiRNAS(gr9)

gr9 <- mobileMethClasses(gr9)

save(gr9,loci9,file="gr9_all.rdata")


round(table(gr9$RdDM[gr9$isTE!="none"])/sum(gr9$isTE!="none"),2)
#
# RdDM_independent        not_known   RdDM_dependent 
#             0.04             0.42             0.54
round(table(gr9$RdDM[gr9$isTE=="none"])/sum(gr9$isTE=="none"),2)

# RdDM_independent        not_known   RdDM_dependent 
#             0.08             0.77             0.15

table(gr9$RdDM[gr9$overlaptype=="SINE"])
# RdDM_independent        not_known   RdDM_dependent 
#                0                7               28

# PolV_independent        not_known   PolV_dependent 
#             0.37             0.29             0.34 
round(table(gr9$polV[gr9$isTE=="none"])/sum(gr9$isTE=="none"),2)

# PolV_independent        not_known   PolV_dependent 
#             0.20             0.71             0.08
