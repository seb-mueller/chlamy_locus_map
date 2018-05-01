#Code written by Nick Matthews to analyse similarities between MCAs with different clusters and dimensions
#Date: 28/01/15

setwd("/home/nem34/segmentMap_II")
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

load("/home/bioinf/nem34/segmentMap_II/ChlamyMCA_Results/MCA100216_4C3D/gr7_clustered.RData")
gr4C3D <- gr7
load("/home/bioinf/nem34/segmentMap_II/ChlamyMCA_Results/MCA110216_5C4D/gr7_clustered.RData")
gr5C4D <- gr7
#load("")
#6C5D <- gr7

comb <- gr4C3D[,1]
comb$gr4C3D <- gr4C3D$cluster
#comb$gr4C4D <- gr4C4D$cluster
comb$gr5C4D <- gr5C4D$cluster
#comb$gr6C4D <- gr6C4D$cluster

#4C3Dvs5C4D
#comb2<-comb[comb$gr4C3Dvs5C4D == "other"]
#new <- rbind(comb2$gr4C3D, comb2$gr5C4D)


comb$gr4C3Dvs5C4D <- rep("other",length(comb))
comb$gr4C3Dvs5C4D[comb$gr4C3D == 1 & comb$gr5C4D == 1] <- "1&1"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 2 & comb$gr5C4D == 1] <- "2&1"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 3 & comb$gr5C4D == 1] <- "3&1"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 4 & comb$gr5C4D == 1] <- "4&1"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 1 & comb$gr5C4D == 2] <- "1&2"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 2 & comb$gr5C4D == 2] <- "2&2"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 3 & comb$gr5C4D == 2] <- "3&2"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 4 & comb$gr5C4D == 2] <- "4&2"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 1 & comb$gr5C4D == 3] <- "1&3"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 2 & comb$gr5C4D == 3] <- "2&3"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 3 & comb$gr5C4D == 3] <- "3&3"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 4 & comb$gr5C4D == 3] <- "4&3"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 1 & comb$gr5C4D == 4] <- "1&4"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 2 & comb$gr5C4D == 4] <- "2&4"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 3 & comb$gr5C4D == 4] <- "3&4"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 4 & comb$gr5C4D == 4] <- "4&4"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 1 & comb$gr5C4D == 5] <- "1&5"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 2 & comb$gr5C4D == 5] <- "2&5"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 3 & comb$gr5C4D == 5] <- "3&5"
comb$gr4C3Dvs5C4D[comb$gr4C3D == 4 & comb$gr5C4D == 5] <- "4&5"

#table(comb$gr4C3Dvs5C4D) 
#  1&1   1&2   1&3   1&5   2&1   2&2   2&3   2&4   3&1   3&2   3&3   3&4   3&5
#11028   233   136     1   232  2273   789     1    16    25    68  1108    13
#  4&1   4&2   4&3   4&4   4&5
#   60    12  1403    17  1685



#table(comb$gr5C4D)
#    1     2     3     4     5
#11336  2543  2396  1126  1699



#Plotting cluster origins
pdf("locusoriginsFinal.pdf")
#par(mfrow=c(2,3))
layout(matrix(c(1,1,2,2,3,3,0,4,4,5,5,0), 2,6,byrow=TRUE))
barplot(c(11336,11028,232,16,60),names.arg=c("Total",1,2,3,4),col=heat.colors(5),main="Cluster 1",ylab="No. of Loci",xlab="Cluster Origin")
barplot(c(2543,233,2273,25,12),names.arg=c("Total",1,2,3,4),col=heat.colors(5),main="Cluster 2",ylab="No. of Loci",xlab="Cluster Origin")
barplot(c(2396,136,789,68,1403),names.arg=c("Total",1,2,3,4),col=heat.colors(5),main="Cluster 3",ylab="No. of Loci",xlab="Cluster Origin")
barplot(c(1126,0,1,1108,17),names.arg=c("Total",1,2,3,4),col=heat.colors(5),main="Cluster 4",ylab="No. of Loci",xlab="Cluster Origin")
barplot(c(1699,1,0,13,1685),names.arg=c("Total",1,2,3,4),col=heat.colors(5),main="Cluster 5",ylab="No. of Loci",xlab="Cluster Origin")
dev.off()

