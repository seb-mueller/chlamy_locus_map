library(FactoMineR)
library(grid)
library(igraph)


#####Setup directories#####
MCAOutputs <- "LociRun2018_multi200_gap100_90c7213_MCAOutputs_05c5bb5"
lociRun <- "LociRun2018_multi200_gap100_90c7213"
baseDir <- "/projects/nick_matthews"
#baseDir <- "C:/Users/Nick/Documents/PhD/Projects/Chlamy"
# Specify location of annotation outputs
inputLocation <- file.path(baseDir, "segmentation_2018", MCAOutputs)
lociLocation <- file.path(baseDir, "segmentation_2018", lociRun)
annoDir     <- file.path(baseDir, "resources")
inputFile <- "gr_fdr0.05.RData"
#Create location for saving figures
figLocation <- file.path(inputLocation,"figures")
try(dir.create(figLocation))
#inputFile <- "gr_fdr0.05.RData"
#set working directory to github repository on cluster
#gitdir      <- file.path(baseDir, "chlamy_locus_map")
gitdir      <- file.path(baseDir, "chlamy_locus_map_github")
#Source code
source(file.path(gitdir,"Scripts/chlamy_source_code.R"))


#Load in files from plotting
load(file.path(inputLocation,"mc6.RData"))
load(file.path(inputLocation,"dimStab.RData"))
load(file.path(inputLocation,"gapStat.RData")); gapStat <- do.call("rbind", gapStat)
load(file.path(inputLocation,"dimList.RData"))
nclust <- 6; ndim <- 7; klist <- dimList[[ndim]]
load(file.path(inputLocation,"clusterings.RData"))
load(file.path(lociLocation,inputFile))


#####Plot large stability plot#####
pdf(file.path(figLocation,"stabilityPlotLarge.pdf"), width = 12, height = 12)
par(mfrow = c(16, 14))
par(mar = c(0.2,0.2,0.2,0.2), oma = c(5,5,5,5))
x0 <- y0 <- Inf; x1 <- y1 <- -Inf
for(ii in 1:16)
    for(jj in 2:15) {
        if(is.matrix(dimStab[[ii]][[jj]])) {
            boxplot((dimStab[[ii]][[jj]]), axes = FALSE, ann = FALSE, ylim = c(0,1))
            box(which = "plot")
        } else plot(NA, NA, xlim = c(0,1), ylim = c(0,1), ylab = "", xlab = "")
        x0 <- min(x0, grconvertX(0, "nfc", "ndc"))
        x1 <- max(x1, grconvertX(1, "nfc", "ndc"))
        y0 <- min(y0, grconvertY(0, "nfc", "ndc"))
        y1 <- max(y1, grconvertY(1, "nfc", "ndc"))
    }
grid.text("# of clusters", x = 0.5, y = y0 / 2)
grid.segments(x0, y0 * 0.75, x1, y0 * 0.75, arrow = arrow(length = unit(0.1, "inches")))
grid.text("Dimensions of MCA used", x = x0 / 2, y = 0.5, rot = 270)
grid.segments(x0 * 0.75, y1, x0 * 0.75, y0, arrow = arrow(length = unit(0.1, "inches")))
dev.off()


#####Variation explained by different dimensions#####
pdf(file.path(figLocation,"VarianceDimensions.pdf"))
par(mar = c(3,4,4,2) + 0.1)
plot(mc6$eig[,2], type = "b", xlab = "Dimension", ylab  = "% of variance", las = 1)
dev.off()

#####Smaller stability plot#####
stab <- c(10, 10)
pdf(file.path(figLocation,"StabilityPlotSmall.pdf"), height = 5, width = 5)
par(mfrow = c(10, 9))
par(mar = c(0.2,0.2,0.2,0.2), oma = c(5,5,5,5))
x0 <- y0 <- Inf; x1 <- y1 <- -Inf
for(ii in 1:stab[1])
    for(jj in 2:stab[2]) {
        if(is.matrix(dimStab[[ii]][[jj]])) {
            boxplot((dimStab[[ii]][[jj]]), axes = FALSE, ann = FALSE, ylim = c(0,1))
            box(which = "plot")
        } else plot(NA, NA, xlim = c(0,1), ylim = c(0,1), ylab = "", xlab = "")
        x0 <- min(x0, grconvertX(0, "nfc", "ndc"))
        x1 <- max(x1, grconvertX(1, "nfc", "ndc"))
        y0 <- min(y0, grconvertY(0, "nfc", "ndc"))
        y1 <- max(y1, grconvertY(1, "nfc", "ndc"))
    }
grid.text("# of clusters", x = 0.5, y = y0 / 2)
grid.segments(x0, y0 * 0.75, x1, y0 * 0.75, arrow = arrow(length = unit(0.1, "inches")))
grid.text("Dimensions of MCA used", x = x0 / 2, y = 0.5, rot = 270)
grid.segments(x0 * 0.75, y1, x0 * 0.75, y0, arrow = arrow(length = unit(0.1, "inches")))
dev.off()

#####Observed vs expected SSE for different cluster numbers#####
pdf(file.path(figLocation,"SumSquares.pdf"), height = 5, width = 5)
plot(x = 2:15, gapStat[,1], ylim = c(6.5,10), ylab = "Wss", xlab = "k", type = "b", col = "red", las = 1)
lines(x = 2:15, gapStat[,2], ylim = c(6.5,10), ylab = "Wss", xlab = "k", type = "b")
legend(x = "bottomleft", lty = 1, col = c("black", "red"), legend = c("Observed", "Expected"), bty = "n")
dev.off()

#####NMI compared to annotation features####
#Need to add some annotations to gr...
annotenv <- new.env()
load(file.path(annoDir,"chlamy_all_annotations.Rdata"), envir = annotenv, verbose = FALSE)
gr <- featureAnn(gr,annotations = annotenv)


pdf(file.path(figLocation,"NMIPlots.pdf"), height = 5, width = 5)
par(mfrow = c(2,1), mar = c(2,3,3,3))
plot(x = 2:length(klist), y = sapply(klist[-1], function(x) igraph::compare(x$cluster, as.integer(as.factor(gr$transposonType)), method = "nmi")), type = "b", xlab = "k", ylab = "NMI", main = "NMI: TE superfamily")

plot(x = 2:length(klist), y = sapply(klist[-1], function(x) igraph::compare(x$cluster, as.integer(as.factor(gr$overlapType)), method = "nmi")), type = "b", xlab = "k", ylab = "NMI", main = "NMI: Annotation features")
dev.off()

pdf(file.path(figLocation,"NMIPlots_NoTE.pdf"), height = 5, width = 5)
plot(x = 2:length(klist), y = sapply(klist[-1], function(x) igraph::compare(x$cluster, as.integer(as.factor(gr$overlapType)), method = "nmi")), type = "b", xlab = "k", ylab = "NMI", main = "NMI: Annotation features")
dev.off()


#####Plot cluster hierarchy#####
z <- c(6,5,4,3,2)

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



