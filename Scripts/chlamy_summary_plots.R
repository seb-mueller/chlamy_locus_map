#Script to run MCA to cluster loci according to their annotations
#Adapted by Nick Matthews from Tom Hardcastle's orginal code
#Date: 13/12/18

#To submit this script to condor
#		 /scripts/conscriptoR /projects/nick_matthews/chlamy_locus_map_github/scripts/chlamy_MCA.r -p19

try(library(FactoMineR))
try(library(clv))
try(library(grid))
library(ggplot2)
try(library(xtable))
try(library(rtracklayer))
try(library(reshape))
library(segmentSeq)
try(library(pROC))
try(library(MASS))
library(RColorBrewer)
try(library(mclust))
try(library(readr))
library(dplyr)

try(library(LabelCompare))
try(library(Kendall))
try(library(RANN))
try(library(igraph))


#####Setup directories#####
MCAOutputs <- "LociRun2018_multi200_gap100_90c7213_MCAOutputs_05c5bb5"
lociRun <- "LociRun2018_multi200_gap100_90c7213"
baseDir <- "/projects/nick_matthews"
#baseDir <- "C:/Users/Nick/Documents/PhD/Projects/Chlamy"
# Specify location of annotation outputs
inputLocation <- file.path(baseDir, "segmentation_2018", MCAOutputs)
lociLocation <- file.path(baseDir, "segmentation_2018", lociRun)
figLocation <- file.path(inputLocation,"figures")
try(dir.create(figLocation))
inputFile <- "gr_fdr0.05.RData"
annoDir     <- file.path(baseDir, "resources")
#set working directory to github repository on cluster
#gitdir      <- file.path(baseDir, "chlamy_locus_map")
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")

#Source code
source(file.path(gitdir,"Scripts/chlamy_source_code.R"))

#Load in loci and gr files
#load(file.path(lociLocation,inputFile))
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05_41c2431.RData")
#Load in list of factors
#factorMaster <- read.csv(file.path(gitdir,"Annotation2Use.csv"),stringsAsFactors = FALSE)

#Specify clusters and dimensions used
nclust <- 6; ndim <- 7

#Load in other files
load(file.path(inputLocation,"clusterings.RData"))
load(file.path(inputLocation,"gr_clustered.RData"))
load(file.path(inputLocation,"resMCA.RData"))

#####Cluster hierarchy#####
#Plots a nice diagram of the hierarchy of clusters
z <- c(6,3,2)

pdf(file.path(figLocation,"cluster_hierarchy.pdf"), height = 7, width = 7)
par(mar = c(1,1,1,1))
plot(NA, NA, ylim = c(0, length(z) * 2 - 1), xlim = c(0, length(clusterings[[5]])), axes = FALSE, xlab = "", ylab = "")
par(mfrow = c(1, 1))
lowTab <- NULL
for(ii in 1:length(z)) {
  
  clid <- z[ii] - 1  
  
  uppTab <- table(clusterings[[clid]])
  uppCoord <- cbind(c(0, cumsum(uppTab)[-length(uppTab)]), cumsum(uppTab))
  
  spacing <- c(0, rep(1 / clid * max(uppCoord) * 0.05, clid))
  
  rect(xleft = uppCoord[,1] + spacing, ybottom = ii * 2 - 2, xright = uppCoord[,2], ytop = ii * 2 - 1, col = rainbow(clid + 1), border = NA)
  #    if(z[ii] == 9)
  #        text(x = (uppCoord[,1] + spacing + uppCoord[,2]) / 2, y = ii * 2 - 1.5, paste0("LC", 1:9))
  
  if(!is.null(lowTab)) {
    normTab <- table(clusterings[[z[ii - 1] - 1]], clusterings[[clid]]) / as.vector(table(clusterings[[z[ii - 1] - 1]]))        
    whCon <- which(normTab > 0.25, arr.ind = TRUE)
    apply(cbind(whCon, normTab[whCon]), 1, function(x)
      igraph:::igraph.Arrows(
        x1 = mean(lowCoord[x[1],]),
        y1 = ii * 2 - 3,
        x2 = mean(uppCoord[x[2],]),
        y2 = ii * 2 - 2, sh.lwd = 10 * x[3]^2, size = 0.2, curve = -0.00005, sh.col = rainbow(z[ii - 1])[x[1]])
    )
  }
  
  lowTab <- uppTab; lowCoord <- uppCoord    
}
dev.off()





#####Chromosome tracks#####
#Load in annotations
load(file.path(annoDir,"chlamy_all_annotations.Rdata"), verbose = FALSE)

#Genes

#mRNA <- genes[genes$type=="mRNA"]
mRNAu <- mRNA[!duplicated(unlist(mRNA$Parent)),]
annottrack_genes <- data.frame(chrom=seqnames(mRNAu), start=start(mRNAu), annot=rep("genes",length(mRNAu)))

#Transposons
annottrack_TEs <- data.frame(chrom=seqnames(transposons), start=start(transposons), annot=rep("TEs",length(transposons)))

#New methylation - combined into one track
methCG=import.gff3(file.path(annoDir,"meth_data/chlamy_CGmeth.gff3"))
methCHH=import.gff3(file.path(annoDir,"meth_data/chlamy_CHHmeth.gff3"))
methCHG=import.gff3(file.path(annoDir,"meth_data/chlamy_CHGmeth.gff3"))
#Ensure sequence levels are matched to the reference
seqlevels(methCG) <- seqlevels(methCHG) <- seqlevels(methCHH) <- seqlevels(gr)
#Create additional object that merged the three methylation files
methAll <- c(methCG,methCHH,methCHG)
annottrack_meth <- data.frame(chrom=seqnames(methAll), start=start(methAll), annot=rep("meth",length(methAll)))

gr<-gr
annottrack_allloci <- data.frame(chrom=as.factor(paste(as.character(gr@seqnames), sep = "")), start=start(gr), annot=rep("loci",length(gr)))
annottrack_cluster <- data.frame(chrom=as.factor(paste(as.character(gr@seqnames), sep = "")), start=start(gr), annot=paste("LC", gr$cluster, sep = ""))

annottrackdf <- rbind(annottrack_genes,annottrack_TEs,annottrack_meth,
                      annottrack_cluster, annottrack_allloci)

annottrackdf$annot <- factor(annottrackdf$annot, levels=c("genes","TEs","meth",paste("LC", as.character(levels(gr$cluster)), sep = ""),"loci"))

###plot cluster tracks - whole genome
gg <- ggplot(annottrackdf) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw()
ggsave(gg,file=file.path(figLocation,"Clustercoverage.eps"),width=30,height=5)
ggsave(gg,file=file.path(figLocation,"Clustercoverage.png"),width=30,height=5)

###plot cluster tracks
#For each of the chromosomes
#Run for all 17 chromosomes...
for(ii in 1:17) {
  annottrackdf_chrX <- subset(annottrackdf, chrom == paste0("chromosome_",ii))
  gg <- ggplot(annottrackdf_chrX) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
  ggsave(gg,file=file.path(figLocation,paste0("Clustercoverage_chr",ii,".eps")),width=10,height=5)
  ggsave(gg,file=file.path(figLocation,paste0("Clustercoverage_chr",ii,".png")),width=10,height=5)
}


#####Paragons#####
#Output paragons (most representative loci for each cluster) for plotting in genome viewer
lapply(1:nclust, function(ii) {
  x <- resMCA$desc.ind$para[[ii]]
  write.table(as.data.frame(gr[as.integer(names(x)),])[,1:4], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA, file = file.path(figLocation,paste("paragons_LC", ii, ".txt", sep = "")))
})

#####Summary tables#####
#Compute summary tables
chrLociTable <- chromosomeLociTable(gr)
write.csv(chrLociTable,file=file.path(figLocation,"chrLociTable.csv"))
annotTable <- annotationSummaryTable(gr,numLevels = 3)
write.csv(annotTable,file=file.path(figLocation,"annotTable.csv"))
clusterAnnotTable <- clusterAnnotationTable(gr)  
write.csv(clusterAnnotTable,file=file.path(figLocation,"clusterAnnotTable.csv"))



