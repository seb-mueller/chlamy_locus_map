

setwd("/home/nem34/segmentMap_II")

# I think most of these libraries aren't necessary anymore. Feel free to try without them and load when you get errors.

#library(clv)
#library(LabelCompare)
library(xtable)
library(rtracklayer)
library(reshape)
library(segmentSeq)
#library(org.At.tair.db)
library(pROC)
library(MASS)
library(RColorBrewer)
#library(Kendall)
#library(RANN)
library(mclust)
library(baySeq)
library(MKmisc)
library(simpleboot)
source("annotation_functions_edit6.r")

# Load up segmentation
load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/segD_first_chlamy_segmentation_nick.RData")

# Select loci based on some fdr
loci7 <- selectLoci(segD, FDR = 0.1, perReplicate = TRUE) #

## creating/exporting coordinates object
gr7 <- loci7@coordinates 

# name the loci 'CR' - chlamydomonas reinhardtii, 'SL' - srna locus
names(gr7) <- sprintf("CRSL%05.f0",1:length(gr7))

#load("gr7_just_cb.RData") #load gr7 you want to add to...

# export as gff3 file for viewing in browser
#export.gff3(gr7,con="loci7_fdr01.gff3")
#import.gff3("loci7_fdr01.gff3")

# the next set of functions are found 'construct_annotation_functions.R'. They take the annotated locus object 'gr7' and add some more annotation data. Each of these functions will need looking at and likely tweaking for chlamy.

# annotate by size class
gr7 <- sizeClasses(gr7)

# annotate with overlapping features
gr7 <- featureAnn(gr7, loci7)

# annotate with counting biases; i.e, is there a higher than average ratio of 21s to 20s, or a higher number of reads starting with As than usual
cl <- makeCluster(24)
gr7 <- countingBiases(gr7,cl)
stopCluster(cl)

#Old methylation function - may still be usefull, picks up a lot more methylation
gr7<-methylation(gr7)

#New methylation functions
gr7 <- methAnnotate(gr7)
gr7 <- methDiff(gr7) #Almost no results - very small datasets

#Extra Current annotations
gr7<-strainSpec(gr7,loci7)
gr7<-lifeCycle(gr7,loci7)
gr7<-mutantSpec(gr7,loci7)
gr7 <- phaseMatch(gr7)


#Other stuff used in Arabidopsis
# this stuff may not be relevant to the chlamy data, so I haven't annotated it. Let me know when you've got through the stuff above, and we can look at where to go next.

#gr7 <- histoneAnnotate(gr7)

#gr7 <- phaseMatch(gr7)

#cl <- makeCluster(24)
#gr7 <- methAnnotate(gr7, cl)
#stopCluster(cl)

#gr7 <- tissueSpec(gr7, loci7) #adapt for zygotes, key mutants, strains
#gr7 <- agoIP(gr7, loci7) #adapt for N starvation


#cl <- makeCluster(24)
#Pol45(gr7, cl)
#stopCluster(cl)

#gr7 <- annPol(gr7)

save(gr7,loci7,file="gr7_all.rdata")

