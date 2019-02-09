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

#Load in loci and gr files
load(file.path(lociLocation,inputFile))
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05_41c2431.RData")
#Load in list of factors
factorMaster <- read.csv(file.path(gitdir,"Annotation2Use.csv"),stringsAsFactors = FALSE)

#####Establish output files#####
#gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
#saveLocation <- file.path(inputLocation, paste(lociRun,"MCAOutputs", gitfingerprint, sep = "_"))
#try(dir.create(saveLocation))


# selected factors which will be used to inform the clustering 
#dplyr and MASS select functions clash so have to specify
selFac <- factorMaster %>% filter(PrimaryAnno==TRUE) %>% dplyr::select(annotation) %>% unlist()

# supplementary factors for which association with clusters will be calculated, but which will not inform the clustering
supFac <- factorMaster %>% filter(SupAnno==TRUE) %>% dplyr::select(annotation) %>% unlist()

#Summary dataframe with the select and supplementary factors
cF6 <- as.data.frame(elementMetadata(gr[,c(selFac,supFac)]))
cF6 <- as.data.frame(unclass(cF6))

# HCPC code from FactoMiner needs tweak to work on kmeans only.
source(file.path(gitdir,"Scripts/hcpc.R"))
#source("C:/Users/Nick/Documents/PhD/Projects/Chlamy/chlamy_locus_map/Scripts/hcpc.R")

#MCA with clusters and dimensions set according to images
#TODO decide cluster and dimension number from the plots
nclust <- 6; ndim <- 7
#nclust <- 2; ndim <- 7
#nclust <- 5; ndim <- 4
#nclust <- 6; ndim <- 4
mc6 <- MCA(cF6, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF6) %in% supFac))
resMCA <- HCPC(mc6, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = nclust, method = "centroid")
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
pvalmatrix[pvalmatrix > 10] <- 10
pvalmatrix[pvalmatrix < -10] <- -10

#Calculates highly correlated pairs
# corp <- matrix(nrow = nrow(pvalmatrix), ncol = nrow(pvalmatrix))
# for(ii in 1:nrow(pvalmatrix)) 
#   for(jj in 1:nrow(pvalmatrix)) corp[ii,jj] <- cor(pvalmatrix[ii,], pvalmatrix[jj,])
# highCor <- which(abs(corp) > 0.9, arr.ind = TRUE)
# highCor <- highCor[!duplicated(t(apply(highCor, 1, sort))),]
# highCor <- highCor[highCor[,1] != highCor[,2],]
# cbind(rownames(pvalmatrix)[highCor[,1]], rownames(pvalmatrix)[highCor[,2]])


colnames(pvalmatrix) <- paste(colnames(pvalmatrix), " (", sapply(colnames(pvalmatrix), function(ii) sum(gr$cluster == ii)), ")", sep = "")


#Select what things we want plotted
selectP <- c(
  paste("sizeclass=", levels(cF6$sizeclass), sep = ""),
  paste("predominant_5prime_letter=", c("A","C","CG","G","T"), sep = ""),
  paste("predominant_sRNA_sizeClass=", unique(cF6$predominant_sRNA_sizeClass),sep=""),
  paste("ratio_strand_class=", levels(cF6$ratio_strand_class), sep = ""),
  paste("phaseClass=phaseClass_", levels(cF6$phaseClass), sep = ""),
  paste("repetitivenessClass=repetitivenessClass_", levels(cF6$repetitivenessClass), sep = ""),
  "methAll=methAll_TRUE",
  "TE_Class_RET=TE_Class_RET_TRUE",
  "TE_Class_DNA=TE_Class_DNA_TRUE",
  "irs=irs_TRUE",
  "trs=trs_TRUE",
  "DCL3dependent=DCL3dependent_TRUE",
  "AGO3dependent=AGO3dependent_TRUE",
  "miRNA=miRNA_TRUE",
  "exons=exons_TRUE","introns=introns_TRUE",
  "promoter=promoter_TRUE",
  "intergenic=intergenic_TRUE",
  paste("expressionClass=", levels(cF6$expressionClass),sep =""),
  "vegetativespecific=vegetativespecific_TRUE",
  "zygotespecific=zygotespecific_TRUE",
  "CC125specific=CC125specific_TRUE",
  "CC1883specific=CC1883specific_TRUE",
  "CC4350specific=CC4350specific_TRUE",
  "Jspecific=Jspecific_TRUE"
)


#rownames(pvalmatrix) <- sapply(strsplit(rownames(pvalmatrixsub),"="),function(x) x[2])

pvalmatrixsel <- pvalmatrix[selectP,]
pvalmatrixsel <- pvalmatrixsel[rowSums(abs(pvalmatrixsel) > 5) > 0,]

#pvalmatrixsel <- pvalmatrixsel[-grep("(^tRNA=)|(^rRNA=)|(^pseudogene=)|(^miRNA=)|(^ncRNA=)|(^lincRNA=)", rownames(pvalmatrixsel)),]
pvalmatrixsub <- pvalmatrix[which(rowSums((pvalmatrix) > 5) > 1),]

#Function to clean up names for plotting
cleanNames <- function(nammat) {    
  #if(length(grep("FALSE", nammat)) > 0) pmat <- pmat[-grep("FALSE", nammat),]
  #if(length(grep("overlaptype=", nammat)) > 0) pmat <- pmat[-grep("overlaptype=", nammat),]
  #if(length(grep("not_known", nammat)) > 0) pmat <- pmat[-grep("not_known", nammat),]
  
  for(ii in 1:length(nammat)) {
    nammat[ii] <- gsub(paste("=", gsub("=.*", "", nammat[ii]), "_", sep = ""), "=", nammat[ii])
    nammat[ii] <- gsub("_new", "", nammat[ii])
    nammat[ii] <- gsub("^has", "", nammat[ii])
    nammat[ii] <- gsub("Class", "", nammat[ii])
    nammat[ii] <- gsub("_class", "", nammat[ii])
    nammat[ii] <- gsub("class", "", nammat[ii])
    nammat[ii] <- gsub("=TRUE", "", nammat[ii])
    nammat[ii] <- gsub("=", ":", nammat[ii])
    nammat[ii] <- gsub("_", " ", nammat[ii])
  }
  nammat <- gsub("median", "moderate", nammat)
  nammat <- gsub("med", "moderate", nammat)
  nammat <- gsub("isIR", "IR", nammat)
  nammat
}

cleanF9 <- cF6; colnames(cleanF9) <- cleanNames(colnames(cF6))
cleanF9 <- cbind(as.data.frame(gr)[,1:3], cleanF9)
#cleanF9 <- cleanF9[,-grep("mobile|mobType|AGO", colnames(cleanF9))]
cleanF9 <- cbind(LC = gr$cluster, cleanF9)

write.table(cleanF9, col.names = NA, quote = FALSE, sep = "\t", file = file.path(inputLocation,"Loci_annotation.txt"))

#Clean up the names
pvalmatrixselC <- pvalmatrixsel
pvalmatrixsubC <- pvalmatrixsub
rownames(pvalmatrixselC) <- cleanNames(rownames(pvalmatrixsel))
rownames(pvalmatrixsubC) <- cleanNames(rownames(pvalmatrixsub))

#Plot heatmap
source(file.path(gitdir,"Scripts/heatmap_centred.R"))
pdf(file.path(figLocation,"featureMatrix_gr6.pdf"),height = 12, width = 10)
heatmap.2(pvalmatrixsubC,
          Rowv = TRUE, Colv = NA, key = FALSE,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(8, 18), cexRow = 1.1, lhei = c(0.01, 5), lwid = c(0.5, 1))
dev.off()


#Plot heatmap of selected features
pdf(file.path(figLocation,"featureMatrixb_gr6.pdf"),width=8, height = 12)
heatmap.2(pvalmatrixselC,
          Colv = NA, Rowv = NA, key = FALSE,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(8, 18), cexRow = 1.1, lhei = c(0.01, 5), lwid = c(0.1, 1))
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

save(gr, file=file.path(inputLocation,"gr_clustered.RData"))
save(resMCA, file=file.path(inputLocation,"resMCA.RData"))
#load("resMCA.RData")

