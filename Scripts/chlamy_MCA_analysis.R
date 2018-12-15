#Script to compute computes diagnostic plots and figures for MCA and clustering
#Adapted by Nick Matthews from Tom Hardcastle's orginal code
#Date: 13/12/18

#To submit this script to condor
#		 /scripts/conscriptoR /projects/nick_matthews/chlamy_locus_map_github/scripts/chlamy_MCA_analysis.r -p19


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
lociRun <- "LociRun2018_multi200_gap100_90c7213"
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
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05.RData")
#load("C:/Users/Nick/Documents/PhD/Projects/Chlamy/gr_fdr0.05_41c2431.RData")
#Load in list of factors
factorMaster <- read_csv(file.path(gitdir,"Annotation2Use.csv"))

#####Establish output files#####
gitfingerprint <- system2("git", args = "rev-parse --short HEAD", stdout = TRUE)
saveLocation <- file.path(inputLocation, paste(lociRun,"MCAOutputs", gitfingerprint, sep = "_"))
try(dir.create(saveLocation))


#####Run Analysis Plots#####
##MCA for categorical data!
#multivariate methods that allows us to analyze the systematic patterns of variations with categorical data
#keep in mind that MCA applies to tables in which the observations are described by a set of qualitative (i.e. categorical) variables

# selected factors which will be used to inform the clustering
#TODO decide exactly which annotations are going in main and supplementary factors

selFac <- factorMaster %>% filter(PrimaryAnno==TRUE) %>% pull(annotation)

# supplementary factors for which association with clusters will be calculated, but which will not inform the clustering
supFac <- factorMaster %>% filter(SupAnno==TRUE) %>% pull(annotation)

#Summary dataframe with the select and supplementary factors
cF7 <- as.data.frame(elementMetadata(gr[,c(selFac,supFac)]))                                                

#png("catchall.png")

#tables summarising output
write("", file = file.path(saveLocation,"ClassTable_gr.txt"))
for(ii in 1:ncol(cF7)) {
    cat(colnames(cF7)[ii], "\t", paste(levels(cF7[,ii]), collapse = "\t"), "\n", file = file.path(saveLocation,"ClassTable_gr.txt"), append = TRUE)
    cat("", "\t", paste(as.numeric(table(cF7[,ii])), collapse = "\t"), "\n\n", file = file.path(saveLocation,"ClassTable_gr.txt"), append = TRUE)
}

# MCA
mc7 <- MCA(cF7, graph = FALSE,ncp = 4 ,quali.sup=which(colnames(cF7) %in% supFac))
#TODO run all the diagnostics
#TODO edit save locations and save names to more appropriate
# parameter sweep on dimensions 1-15 and clusters 2-15
cl <- makeCluster(14)
dimList <- list()
for(nn in 1:15) {
    mc7 <- MCA(cF7, graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF7) %in% supFac))
    dimList[[nn]] <- c(list(NA), parLapply(cl, 2:15, function(i, coords)
        kmeans(x = coords, iter.max = 1000, nstart = 1000, centers = i), coords = mc7$ind$coord))
}
    
save(dimList, file = file.path(saveLocation,"dimList.RData"))
#load("dimList.RData")

# probably not useful as a measure of performance
png(file.path(saveLocation,"kmean_ss.png"))
image(sapply(dimList, function(x) sapply(x[-1], function(y) y$betweenss / y$totss)))
dev.off()

# rand test to compare overlap with transposable element superfamilies 
zzz <- (sapply(dimList, function(x) sapply(x[-1], function(y) adjustedRandIndex(y$cluster, gr$TE))))
png(file.path(saveLocation,"kmean_randTE.png"))
image(zzz)
dev.off()

# stability analyses. Takes a while!
dimStab <- list()
for(nn in 1:15) {
    dimStab[[nn]] <- list()
    for(nclust in 2:15) {
        kcb <- list()
        message(nn, ":", nclust, appendLF = FALSE)
        for(ii in 1:10) {
            repeat {
                message(".", appendLF = FALSE)
                rsamp <- sample(1:nrow(cF7), nrow(cF7), replace = TRUE)
                mcb <- try(MCA(cF7[rsamp,], graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF7) %in% supFac)))
                if(!"try-error" %in% class(mcb)) break()
                
            }
            cob <- mcb$ind$coord
            kcb[[ii]] <- list(rsamp = rsamp,
                        km = kmeans(cob, centers = nclust, iter.max = 1000, nstart = 100),
                        cob = mcb$ind$coord, eig = mcb$eig)
        }

        kjaq <- sapply(kcb, function(x, kc) {
            bov <- (table(cbind.data.frame(kc = kc$cluster[x$rsamp], boot = x$km$cluster)[!duplicated(x$rsamp),]))
            apply(bov / (outer(rowSums(bov) , colSums(bov), FUN='+') - bov), 2, max)
        }, kc = dimList[[nn]][[nclust]])
        dimStab[[nn]][[nclust]] <- kjaq
        message()
    }
}

save(dimStab, file = file.path(saveLocation,"dimStab.RData"))
#load("dimStab.RData")

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

# Davies-Bouldin
dviM <- apply(sapply(dimList, function(x, mc) parLapply(cl, x[-1], function(y, mc)
    clv.Davies.Bouldin(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "complete", "complete")
                                                      , mc = mc), mc = mc7),2,unlist)

dviM2 <- apply(sapply(dimList, function(x, mc)
    zz <- parLapply(cl, x[-1], function(y, mc)
        clv.Davies.Bouldin(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "average", "average")
                  , mc = mc), mc = mc7),2,unlist)

dviM3 <- apply(sapply(dimList, function(x, mc)
    zz <- parLapply(cl, x[-1], function(y, mc)
        clv.Davies.Bouldin(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "centroid", "centroid")
                  , mc = mc), mc = mc7),2,unlist)

pdf(file.path(saveLocation,"davies_bouldin_image.pdf"))
image(dviM)
image(dviM2)
image(dviM3)
dev.off()

# Dunn
dunnM <- apply(sapply(dimList, function(x, mc) parLapply(cl, x[-1], function(y, mc)
    clv.Dunn(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "complete", "complete")
                                                      , mc = mc), mc = mc7),2,unlist)

dunnM2 <- apply(sapply(dimList, function(x, mc)
    zz <- parLapply(cl, x[-1], function(y, mc)
        clv.Dunn(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "average", "average")
                  , mc = mc), mc = mc7),2,unlist)

dunnM3 <- apply(sapply(dimList, function(x, mc)
    zz <- parLapply(cl, x[-1], function(y, mc)
        clv.Dunn(cls.scatt.data(mc$ind$coord, as.integer(y$cluster)), "centroid", "centroid")
                  , mc = mc), mc = mc7),2,unlist)

pdf(file.path(saveLocation,"dunn_image.pdf"))
image(dunnM)
image(dunnM2)
image(dunnM3)
dev.off()

#Examine images and choose correct dimension/cluster number. Need to talk to Tom for interpretation
