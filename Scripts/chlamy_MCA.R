#Script to run MCA to cluster loci according to their annotations
#Adapted by Nick Matthews from Tom Hardcastle's orginal code
#Date: 13/12/18

#To submit this script to condor
#		 /scripts/conscriptoR /projects/nick_matthews/chlamy_locus_map_github/scripts/chlamy_MCA.r -p19

library(FactoMineR)
library(clv)
library(grid)
library(ggplot2)
library(xtable)
library(rtracklayer)
library(reshape)
library(segmentSeq)
library(pROC)
library(MASS)
library(RColorBrewer)
library(mclust) 
library(readr)
library(dplyr)

#####Setup directories#####
lociRun <- "LociRun2018_multi200_gap100"
baseDir <- "/projects/nick_matthews"
#baseDir <- "C:/Users/Nick/Documents/PhD/Projects/Chlamy"
# Specify location of annotation outputs
inputLocation <- file.path(baseDir, "segmentation_2018", lociRun)
inputfile <- "gr_fdr0.05.RData"
#set working directory to github repository on cluster
#gitdir      <- file.path(baseDir, "chlamy_locus_map")
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")

#Load in loci and gr files
load(file.path(inputLocation,inputFile))
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05_41c2431.RData")
#Load in list of factors
factorMaster <- read_csv(file.path(gitdir,"Annotation2Use.csv"))

#####Establish output files#####
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
saveLocation <- file.path(inputLocation, paste(lociRun,"MCAOutputs", gitfingerprint, sep = "_"))
try(dir.create(saveLocation))


# selected factors which will be used to inform the clustering
selFac <- factorMaster %>% filter(PrimaryAnno==TRUE) %>% pull(annotation)

# supplementary factors for which association with clusters will be calculated, but which will not inform the clustering
supFac <- factorMaster %>% filter(SupAnno==TRUE) %>% pull(annotation)

#Summary dataframe with the select and supplementary factors
cF7 <- as.data.frame(elementMetadata(gr[,c(selFac,supFac)]))    

# HCPC code from FactoMiner needs tweak to work on kmeans only.
source(file.path(gitdir,"Scripts/hcpc.R"))
#source("C:/Users/Nick/Documents/PhD/Projects/Chlamy/chlamy_locus_map/Scripts/hcpc.R")

#MCA with clusters and dimensions set according to images
#TODO decide cluster and dimension number from the plots
nclust <- 4; ndim <- 3
#nclust <- 4; ndim <- 4
#nclust <- 5; ndim <- 4
#nclust <- 6; ndim <- 4
mc7 <- MCA(cF7, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF7) %in% supFac))
resMCA <- HCPC(mc7, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = nclust, method = "centroid")
gr$cluster <- as.factor(resMCA$data.clust$clust)

# now check feature associations


# extract names
rownames(resMCA$call$t$res$var$v.test)
selectionCat <- unique(do.call("c", lapply(resMCA$desc.var$category, row.names)))
selectionCat <- selectionCat[-grep("\\.NA",selectionCat)]


# extract v-test values
vtestmatrix <- matrix(nr=length(selectionCat),ncol=nclust)
rownames(vtestmatrix) <- selectionCat
colnames(vtestmatrix) <- paste("",1:nclust,col="",sep="")

# extract p-values
pvalmatrix <- vtestmatrix
for (i in 1:nclust) {
  tmp <- resMCA$desc.var$category[[i]]
  # rownames(tmp) <- sapply(strsplit(rownames(tmp),"="),function(x) x[2])
  vtestmatrix[,i] <- log2(tmp[selectionCat,"Mod/Cla"]/tmp[selectionCat,"Global"])
  pvalmatrix[,i] <- tmp[selectionCat,"p.value"]
}

# map to log-scale and set dynamic range
pvalmatrix <- log10(pvalmatrix)
pvalmatrix[vtestmatrix < 0] <- -pvalmatrix[vtestmatrix < 0]
pvalmatrix[pvalmatrix > 255] <- 255
pvalmatrix[pvalmatrix < -255] <- -255


colnames(pvalmatrix) <- paste(colnames(pvalmatrix), " (", sapply(colnames(pvalmatrix), function(ii) sum(gr$cluster == ii)), ")", sep = "")

# select only significant classes
pvalmatrixsub <- pvalmatrix[which(rowSums(abs(pvalmatrix) > 3) > 1),]
pvalmatrixsub <- pvalmatrixsub[-grep("FALSE", rownames(pvalmatrixsub)),]
#pvalmatrixsub <- pvalmatrixsub[-grep("overlaptype=", rownames(pvalmatrixsub)),]

#TODO decide the select variables to use, and modify plot to have nice row names
# or select specific named classes
selectP <- c(
  paste("sizeclass=", levels(cF7$sizeclass), sep = ""),
  "ratio21vs20Class=ratio21vs20Class_low","ratio21vs20Class=ratio21vs20Class_med","ratio21vs20Class=ratio21vs20Class_high",
  "ratioSmallvsNormalClass=ratioSmallvsNormalClass_low","ratioSmallvsNormalClass=ratioSmallvsNormalClass_med","ratioSmallvsNormalClass=ratioSmallvsNormalClass_high",
  "ratioBigvsNormalClass=ratioBigvsNormalClass_low","ratioBigvsNormalClass=ratioBigvsNormalClass_med","ratioBigvsNormalClass=ratioBigvsNormalClass_high",
  paste("ratio_strand_class=", levels(cF7$ratio_strand_class), sep = ""),
  paste("predominant_5prime_letter=", c("A","AC","C","CG","CT","G","GT","T"), sep = ""),
  paste("repetitivenessClass=repetitivenessClass_", levels(cF7$repetitivenessClass), sep = ""),
  paste("expressionClass=", levels(cF7$expressionClass),sep =""),
  paste("Phased=Phased_", levels(cF7$Phased),sep =""),
  "meth=meth_TRUE",
  "methCG=methCG_TRUE",
  "vegetativespecific=vegetativespecific_TRUE",
  "zygotespecific=zygotespecific_TRUE",
  "DCL3dependent=DCL3dependent_TRUE",
  "gene=gene_TRUE",
  "CDS=CDS_TRUE","exon=exon_TRUE","fiveprimeUTR=fiveprimeUTR_TRUE","threeprimeUTR=threeprimeUTR_TRUE","promoter=promoter_TRUE","Jspecific=Jspecific_TRUE",
  "miRNA=miRNA_TRUE",
  "IR=IR_TRUE",
  "TR=TR_TRUE",
  paste("TEclass=TEclass_", c("none","DNA","RET") , sep = ""),
  paste("TEorder=TEorder_", c("LTR","LINE","SINE","TIR") , sep = "")
)


pvalmatrixsel <- pvalmatrix[selectP,]

# heatmap function needs some tweaks to work on these data
source(file.path(gitdir,"Scripts/heatmap_centred.R"))
#source("C:/Users/Nick/Documents/PhD/Projects/Chlamy/chlamy_locus_map/Scripts/heatmap_centred.R")

#examine heatmaps
for(ii in 1:nrow(pvalmatrixsub))
  rownames(pvalmatrixsub)[ii] <- gsub(paste("=", gsub("=.*", "", rownames(pvalmatrixsub)[ii]), "_", sep = ""), "=", rownames(pvalmatrixsub)[ii])
png("featureMatrix4_new_gr.png",width=1200, height = 1600)
par(mar = c(5.1, 4.1, 9.1, 2.1))
heatmap.2(pvalmatrixsub,
          Colv = NA,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(8, 20), cexRow = 1.1)
dev.off()

for(ii in 1:nrow(pvalmatrixsel))
  rownames(pvalmatrixsel)[ii] <- gsub(paste("=", gsub("=.*", "", rownames(pvalmatrixsel)[ii]), "_", sep = ""), "=", rownames(pvalmatrixsel)[ii])
png("featureMatrix4b_gr.png",width=1200, height = 1600)
heatmap.2(pvalmatrixsel,
          Colv = NA, Rowv = NA,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(8, 14), cexRow = 1.1)
dev.off()


# remove boring stuff for write to tables
categories <- resMCA$desc.var$category
filcat <- lapply(categories, function(x) {
  x <- x[x[,5] > 0,]
  if(length(grep("overlaptype=", rownames(x))) > 0) x <- x[-grep("overlaptype=", rownames(x)),]
  if(length(grep("=.*NA", rownames(x))) > 0) x <- x[-grep("=.*NA", rownames(x)),]
  if(length(grep("=.*FALSE", rownames(x))) > 0) x <- x[-grep("=.*FALSE", rownames(x)),]
  if(length(grep("=.*not_known", rownames(x))) > 0) x <- x[-grep("=.*not_known", rownames(x)),]
  if(length(grep("=.*none", rownames(x))) > 0) x <- x[-grep("=.*none", rownames(x)),]
  x
})

save(gr, file=file.path(saveLocation,"gr_clustered.RData"))
save(resMCA, file=file.path(saveLocation,"resMCA.RData"))
#load("resMCA.RData")


# try plotting density of loci in different clusters across genome; compare with gene and transposable element densities.
#TODO decide what's going on the chromosome tracks
#Genes
genes <- import.gff3("/home/bioinf/nem34/Creinhardtii_281_v5.5.gene_exons.gff3")
mRNA <- genes[genes$type=="mRNA"]
mRNAu <- mRNA[!duplicated(unlist(mRNA$Parent)),]
annottrack_genes <- data.frame(chrom=as.factor(as.character(mRNAu@seqnames)), start=start(mRNAu), annot=rep("genes",length(mRNAu)))

#Trnsposons
transposonsgff <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/transposons_newcut.gff3") #GRanges
annottrack_TEs <- data.frame(chrom=as.factor(as.character(transposonsgff@seqnames)), start=start(transposonsgff), annot=rep("TEs",length(transposonsgff)))

#Methylation
meth<-import.gff3("/data/pipeline/prod/SL55/SL55.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff3")
annottrack_meth <- data.frame(chrom=as.factor(as.character(meth@seqnames)), start=start(meth), annot=rep("meth",length(meth)))

#New methylation - combined into one track
methCG<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CGmeth.gff3")
annottrack_methCG <- data.frame(chrom=as.factor(as.character(methCG@seqnames)), start=start(methCG), annot=rep("meth2",length(methCG)))
methCHH<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHHmeth.gff3")
annottrack_methCHH <- data.frame(chrom=as.factor(as.character(methCHH@seqnames)), start=start(methCHH), annot=rep("meth2",length(methCHH)))
methCHG<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHGmeth.gff3")
annottrack_methCHG <- data.frame(chrom=as.factor(as.character(methCHG@seqnames)), start=start(methCHG), annot=rep("meth2",length(methCHG)))
annottrack_meth2 <- rbind(annottrack_methCG,annottrack_methCHH,annottrack_methCHG)

#loci
gr<-gr;prefix="loci_"
annottrack_allloci <- data.frame(chrom=as.factor(paste(as.character(gr@seqnames), sep = "")), start=start(gr), annot=rep("all loci",length(gr)))
annottrack_cluster <- data.frame(chrom=as.factor(paste(as.character(gr@seqnames), sep = "")), start=start(gr), annot=paste("LC", gr$cluster, sep = ""))

annottrackdf <- rbind(annottrack_genes,annottrack_TEs,annottrack_meth,annottrack_meth2,
                      annottrack_cluster, annottrack_allloci)

annottrackdf$annot <- factor(annottrackdf$annot, levels=c("genes","TEs","meth","meth2",paste("LC", as.character(levels(gr$cluster)), sep = ""),"all loci"))

###plot cluster tracks - whole genome
gg <- ggplot(annottrackdf) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw()
ggsave(gg,file=paste0(prefix,"Clustercoverage.eps"),width=30,height=5)
ggsave(gg,file=paste0(prefix,"Clustercoverage.png"),width=30,height=5)

###plot cluster tracks
annottrackdf_chr1 <- subset(annottrackdf, chrom == "chromosome_1")
gg <- ggplot(annottrackdf_chr1) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
ggsave(gg,file="Clustercoverage_chr1.eps",width=10,height=5)
ggsave(gg,file="Clustercoverage_chr1.png",width=10,height=5)

annottrackdf_chr2 <- subset(annottrackdf, chrom == "chromosome_2")
gg <- ggplot(annottrackdf_chr2) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
ggsave(gg,file="Clustercoverage_chr2.eps",width=10,height=5)
ggsave(gg,file="Clustercoverage_chr2.png",width=10,height=5)

annottrackdf_chr3 <- subset(annottrackdf, chrom == "chromosome_3")
gg <- ggplot(annottrackdf_chr3) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
ggsave(gg,file="Clustercoverage_chr3.eps",width=10,height=5)
ggsave(gg,file="Clustercoverage_chr3.png",width=10,height=5)

annottrackdf_chr4 <- subset(annottrackdf, chrom == "chromosome_4")
gg <- ggplot(annottrackdf_chr4) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
ggsave(gg,file="Clustercoverage_chr4.eps",width=10,height=5)
ggsave(gg,file="Clustercoverage_chr4.png",width=10,height=5)

#annottrackdf_chr10 <- subset(annottrackdf, chrom == "chromosome_10")
#gg <- ggplot(annottrackdf_chr10) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
#ggsave(gg,file="Clustercoverage_chr10.eps",width=10,height=5)
#ggsave(gg,file="Clustercoverage_chr10.png",width=10,height=5)

#annottrackdf_sfld18 <- subset(annottrackdf, chrom == "scaffold_18")
#gg <- ggplot(annottrackdf_csfld18) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw(base_size = 8)
#ggsave(gg,file="Clustercoverage_sfld18.eps",width=10,height=5)
#ggsave(gg,file="Clustercoverage_sfld18.png",width=10,height=5)


#Output paragons (most representative loci for each cluster) for plotting in genome viewer
lapply(1:nclust, function(ii) {
  x <- resMCA$desc.ind$para[[ii]]
  write.table(as.data.frame(gr[as.integer(names(x)),])[,1:4], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA, file = paste("paragons_LC", ii, ".txt", sep = ""))
})

dev.off() 