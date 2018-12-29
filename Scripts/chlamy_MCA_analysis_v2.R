#Script to compute computes diagnostic plots and figures for MCA and clustering
#Adapted by Nick Matthews from Tom Hardcastle's orginal code
#Date: 13/12/18

#To submit this script to condor
#		 /scripts/conscriptoR /projects/nick_matthews/chlamy_locus_map_github/scripts/chlamy_MCA_analysis.r -p19


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

library(LabelCompare)
try(library(Kendall))
try(library(RANN))
try(library(igraph))


#####Setup directories#####
lociRun <- "LociRun2018_multi200_gap100_90c7213"
baseDir <- "/projects/nick_matthews"
#baseDir <- "C:/Users/Nick/Documents/PhD/Projects/Chlamy"
# Specify location of annotation outputs
inputLocation <- file.path(baseDir, "segmentation_2018", lociRun)
inputFile <- "gr_fdr0.05.RData"
#set working directory to github repository on cluster
#gitdir      <- file.path(baseDir, "chlamy_locus_map")
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")

#Load in gr file
load(file.path(inputLocation,inputFile))
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05.RData")
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05_41c2431.RData")
#Load in list of factors
factorMaster <- read.csv(file.path(gitdir,"Annotation2Use.csv"),stringsAsFactors = FALSE)

#####Establish output files#####
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
saveLocation <- file.path(baseDir,"segmentation_2018", paste(lociRun,"MCAOutputs", gitfingerprint, sep = "_"))
try(dir.create(saveLocation))


#####Run Analysis Plots#####
##MCA for categorical data!
#multivariate methods that allows us to analyze the systematic patterns of variations with categorical data
#keep in mind that MCA applies to tables in which the observations are described by a set of qualitative (i.e. categorical) variables

# selected factors which will be used to inform the clustering
#TODO decide exactly which annotations are going in main and supplementary factors

selFac <- factorMaster %>% filter(PrimaryAnno==TRUE) %>% select(annotation) %>% unlist()

# supplementary factors for which association with clusters will be calculated, but which will not inform the clustering
supFac <- factorMaster %>% filter(SupAnno==TRUE) %>% select(annotation) %>% unlist()

#Summary dataframe with the select and supplementary factors
cF9 <- as.data.frame(elementMetadata(gr[,c(selFac,supFac)]))                                                
cF9 <- as.data.frame(unclass(cF9))  
#png("catchall.png")

#tables summarising output
write("", file = file.path(saveLocation,"ClassTable_gr.txt"))
for(ii in 1:ncol(cF9)) {
    cat(colnames(cF9)[ii], "\t", paste(levels(cF9[,ii]), collapse = "\t"), "\n", file = file.path(saveLocation,"ClassTable_gr.txt"), append = TRUE)
    cat("", "\t", paste(as.numeric(table(cF9[,ii])), collapse = "\t"), "\n\n", file = file.path(saveLocation,"ClassTable_gr.txt"), append = TRUE)
}

# MCA
mc9 <- MCA(cF9, graph = FALSE,ncp = 4 ,quali.sup=which(colnames(cF9) %in% supFac))
#TODO run all the diagnostics
#TODO edit save locations and save names to more appropriate
# parameter sweep on dimensions 1-15 and clusters 2-15
cl <- makeCluster(14)
dimList <- list()
for(nn in 1:15) {
    mc9 <- MCA(cF9, graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF9) %in% supFac))
    dimList[[nn]] <- c(list(NA), parLapply(cl, 2:15, function(i, coords)
        kmeans(x = coords, iter.max = 1000, nstart = 1000, centers = i), coords = mc9$ind$coord))
}
    
save(dimList, file = file.path(saveLocation,"dimList.RData"))
#load("dimList.RData")
stopCluster(cl)



# stability analyses. Takes a while!
cl <- makeCluster(40)
clusterEvalQ(cl, library(FactoMineR))
dimStab <- list()
for(nn in 1:16) {
  dimStab[[nn]] <- list()
  for(nclust in 2:15) {
    message(nn, ":", nclust, appendLF = FALSE)
    dimStab[[nn]][[nclust]] <- do.call("rbind", parLapply(cl, 1:100, function(ii, kc, nn, nclust, cF9, supFac) {
      message(".", appendLF = FALSE)
      repeat {
        
        rsamp <- unique(sample(1:nrow(cF9), nrow(cF9), replace = TRUE))
        mcb <- try(MCA(cF9[unique(rsamp),], graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF9) %in% supFac)))
        if(!"try-error" %in% class(mcb)) break
      }
      cob <- mcb$ind$coord
      #kcb[[ii]] <- list(rsamp = rsamp,
      km = kmeans(cob, centers = nclust, iter.max = 1000, nstart = 100)
      #            cob = mcb$ind$coord, eig = mcb$eig)
      bov <- (table(cbind.data.frame(kc = kc$cluster[rsamp], boot = km$cluster)))#[!duplicated(rsamp),]))
      mstat <- apply(bov / (outer(rowSums(bov) , colSums(bov), FUN='+') - bov), 2, max)
      return(mstat)
    }, kc = dimList[[nn]][[nclust]], nn = nn, nclust = nclust, cF9 = cF9, supFac = supFac)
    )
    message()
  }
}

save(dimStab, file = file.path(saveLocation,"dimStab.RData"))
#load("dimStab.RData")

stopCluster(cl)

#Do clusterings with identified settings
nclust <- 6; ndim <- 7
klist <- dimList[[ndim]]
mc9 <- MCA(cF9, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF9) %in% supFac))
save(mc9, file = "mc9.RData")


#Calculate gapstat based on Tibshirani et al. 2001 and using Hardcastle et al. 2018 code
cl <- makeCluster(10)
gapStat <- lapply(2:15, function(cls) {
  uW <- parSapply(cl, 1:10, function(rrr, cls, klist, mc9) {
    clsplit <-split(1:nrow(mc9$ind$coord), klist[[cls]]$cluster)
    bbox <- lapply(clsplit, function(z) apply(mc9$ind$coord[z,,drop = FALSE], 2, range))
    rdat <- do.call("rbind", lapply(1:cls, function(clust)
      apply(bbox[[clust]], 2, function(x) runif(sum(klist[[cls]]$cluster == cls), min = x[1], max = x[2]))
    ))
    kmr <- kmeans(x = rdat, iter.max = 1000, nstart = 1000, centers = cls)
    log(sum(kmr$withinss))
  }, cls = cls, klist = klist, mc9 = mc9)
  se <- sd(uW) * sqrt(1 + 1 / length(uW))
  c(mean(uW), log(sum(klist[[cls]]$withinss)), se)
})

save(gapStat, file = file.path(saveLocation,"gapStat.RData"))
stopCluster(cl)

clusterings <- lapply(2:nclust, function(kk) {    
  mc9 <- MCA(cF9, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF9) %in% supFac))
  resMCA <- HCPC(mc9, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = kk, method = "centroid")
  as.factor(resMCA$data.clust$clust)
})
save(clusterings, file = file.path(saveLocation,"clusterings.RData"))



#Plot stability analysis

pdf(file.path(saveLocation,"stabilityplots.pdf"),20,20)
par(mfrow = c(15,15), mar = c(0.2,0.2,0.2,0.2))
for(ii in 1:15)
    for(jj in 1:15)
        if(is.matrix(dimStab[[ii]][[jj]])) {
            boxplot(t(dimStab[[ii]][[jj]]), axes = FALSE, ylim = c(0,1))
        } else plot(NA, NA, xlim = c(0,1), ylim = c(0,1), axes = FALSE, ylab = "", xlab = "")
dev.off()
		
		
# set of tests on internal clustering performance
cl <- makeCluster(19)
clusterEvalQ(cl, library(clv))



#Examine images and choose correct dimension/cluster number. Need to talk to Tom for interpretation
