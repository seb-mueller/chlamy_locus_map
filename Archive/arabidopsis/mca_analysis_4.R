library(FactoMineR)
library(clv)
library(LabelCompare)
library(grid)
library(ggplot2)
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
library(igraph)

source("hcpc.R")
source("/home/sm934/code/R/seb_functions.r")
load("lociIV_fdr01_c.rdata")
load("gr9_all.rdata")

source("~/Code/extractSequences.R")
tair <- extractSequences("/data/public_data/arabidopsis/TAIR_10/TAIR10_chr_all.fa")

load("/home/tjh48/Code/root_grafts_methylation/newFigures/mobileTypes.RData")
names(mobileTypes)
for(mm in 1:length(mobileTypes)) {
    seqlevels(mobileTypes[[mm]]@coordinates) <- as.character(1:5)
    values(gr9)[,paste("mobType_", gsub(">|=|:", ".", names(mobileTypes)[mm]), sep = "")] <- segmentSeq:::.getOverlaps(gr9, mobileTypes[[mm]]@coordinates, whichOverlaps = FALSE)
}

##MCA for categoriacal data!
#multivariate methods that allows us to analyze the systematic patterns of variations with categorical data
#keep in mind that MCA applies to tables in which the observations are described by a set of qualitative (i.e. categorical) variables
selFac <- c("sizeclass", "polIV_new", "h3Class", "polV_new", "rdr2_new", "ratio_strand_class","predominant_5prime_letter","hasExpression",
            "ratio21vs24Class","tissue_expressionClass", "repetitivenessClass","isPhased", "dclTripleClass",
            "hasH3K27me1","hasH3K27me3","hasH3K9me2","hasH3K9Ac","hasH3K4me2","hasH3K4me3","hasH3K36me2",
            "CGclass", "CHGclass", "CHHclass")                       

z <- (values(gr9)[,colnames(values(gr9))[grep("\\.mobile[A-F]", colnames(values(gr9)))]])
mobileColumns <- grep("mobile", colnames(values(gr9)))
supFac <- c("overlaptype","isIR", 
            colnames(values(gr9))[which(colnames(values(gr9)) %in% gsub("/|\\?", ".", levels(gr9$overlaptype)))],
            colnames(values(gr9))[grep("mobType", colnames(values(gr9)))],
            colnames(values(gr9))[grep("\\.mobile[A-F]", colnames(values(gr9)))][sapply(1:18, function(x) sum(as.integer(z[,x]))) > 0],
            "Ago1IPMclass","Ago2IPMclass","Ago4IPMclass","Ago6IPclass","Ago5IPclass","Ago9IPclass")

cF9 <- as.data.frame(elementMetadata(gr9[,c(selFac,supFac)]))                                                      
                                                      
write("", file = "ClassTable_gr9.txt")
for(ii in 1:ncol(cF9)) {
    cat(colnames(cF9)[ii], "\t", paste(levels(cF9[,ii]), collapse = "\t"), "\n", file = "ClassTable_gr9.txt", append = TRUE)
    cat("", "\t", paste(as.numeric(table(cF9[,ii])), collapse = "\t"), "\n\n", file = "ClassTable_gr9.txt", append = TRUE)
}

mc9 <- MCA(cF9, graph = FALSE,ncp = 4 ,quali.sup=which(colnames(cF9) %in% supFac))

cl <- makeCluster(14)
dimList <- list()
for(nn in 1:16) {
    mc9 <- MCA(cF9, graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF9) %in% supFac))
    dimList[[nn]] <- c(list(NA), parLapply(cl, 2:15, function(i, coords)
        kmeans(x = coords, iter.max = 1000, nstart = 1000, centers = i), coords = mc9$ind$coord))
}
save(dimList, file = "dimList_9.RData")

stopCluster(cl)
load("dimList_9.RData")

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

save(dimStab, file = "dimStab_9_noIP.RData")


nclust <- 9; ndim <- 6
klist <- dimList[[ndim]]
mc9 <- MCA(cF9, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF9) %in% supFac))
save(mc9, file = "mc9.RData")

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

save(gapStat, file = "gapStat_dim6.RData")

clusterings <- lapply(2:nclust, function(kk) {    
    mc9 <- MCA(cF9, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF9) %in% supFac))
    resMCA <- HCPC(mc9, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = kk, method = "centroid")
    as.factor(resMCA$data.clust$clust)
})
save(clusterings, file = "clusterings_noIP.RData")

mc9 <- MCA(cF9, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF9) %in% supFac))
resMCA <- HCPC(mc9, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = nclust, method = "centroid")
gr9$cluster <- as.factor(resMCA3$data.clust$clust)
save(gr9, resMCA, supFac, selFac, file = "gr9_clustered_update.rdata")

which(colnames(values(gr9)) %in% c(supFac, selFac))
annotatedLoci <- cbind(as.data.frame(gr9)[,1:4], cluster = values(gr9)$cluster, values(gr9)[,which(colnames(values(gr9)) %in% c(supFac, selFac))])
colnames(annotatedLoci)[1] <- "chr"
write.table(annotatedLoci, file = "~/My_Papers/annotating-plant-epigenome/annotated_sRNA_loci.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# category descriptions for the locus classes
lapply(resMCA$desc.var$category, function(x) x[-grep("FALSE|NA|none|not_known", rownames(x)),])[[3]]    

selectionCat <- unique(do.call("c", lapply(resMCA$desc.var$category, row.names)))
selectionCat <- selectionCat[-grep("\\.NA",selectionCat)]


vtestmatrix <- matrix(nr=length(selectionCat),ncol=nclust)
rownames(vtestmatrix) <- selectionCat
colnames(vtestmatrix) <- paste("",1:nclust,col="",sep="")
pvalmatrix <- vtestmatrix

for (i in 1:nclust) {
    tmp <- resMCA$desc.var$category[[i]]
                                        # rownames(tmp) <- sapply(strsplit(rownames(tmp),"="),function(x) x[2])
    vtestmatrix[,i] <- log2(tmp[selectionCat,"Mod/Cla"]/tmp[selectionCat,"Global"])
    pvalmatrix[,i] <- tmp[selectionCat,"p.value"]
}

pvalmatrix <- log10(pvalmatrix)
pvalmatrix[vtestmatrix < 0] <- -pvalmatrix[vtestmatrix < 0]
pvalmatrix[pvalmatrix > 10] <- 10
pvalmatrix[pvalmatrix < -10] <- -10

corp <- matrix(nrow = nrow(pvalmatrix), ncol = nrow(pvalmatrix))
for(ii in 1:nrow(pvalmatrix)) 
    for(jj in 1:nrow(pvalmatrix)) corp[ii,jj] <- cor(pvalmatrix[ii,], pvalmatrix[jj,])
highCor <- which(abs(corp) > 0.9, arr.ind = TRUE)
highCor <- highCor[!duplicated(t(apply(highCor, 1, sort))),]
highCor <- highCor[highCor[,1] != highCor[,2],]
cbind(rownames(pvalmatrix)[highCor[,1]], rownames(pvalmatrix)[highCor[,2]])


colnames(pvalmatrix) <- paste(colnames(pvalmatrix), " (", sapply(colnames(pvalmatrix), function(ii) sum(gr9$cluster == ii)), ")", sep = "")

h3ids <- colnames(cF9)[grep("hasH3", colnames(cF9))]

selectP <- c(
    paste("sizeclass=", levels(cF9$sizeclass), sep = ""),
    "ratio21vs24Class=ratio21vs24Class_low","ratio21vs24Class=ratio21vs24Class_med","ratio21vs24Class=ratio21vs24Class_high",
    paste("ratio_strand_class=", levels(cF9$ratio_strand_class), sep = ""),
    paste("predominant_5prime_letter=", c("A", "AT", "CGT"), sep = ""),
    paste("polIV_new=", c("polIV_new_PolIV_dependent", "polIV_new_PolIV_independent"), sep = ""),    
    paste("polV_new=", c("polV_new_PolV_dependent", "polV_new_PolV_independent"), sep = ""),
    paste("rdr2_new=", c("rdr2_new_Rdr2_dependent", "rdr2_new_Rdr2_independent"), sep = ""),
    paste("repetitivenessClass=", c("repetitivenessClass_high", "repetitivenessClass_med", "repetitivenessClass_low"), sep = ""),
    paste("tissue_expressionClass=", levels(cF9$tissue_expressionClass),sep =""),
    "isPhased=isPhased_moderate",
    "isIR=isIR_TRUE",
    paste("CGclass=CGclass", levels(cF9$CGclass), sep = "_"),
    paste("CHGclass=CHGclass", levels(cF9$CHGclass), sep = "_"),
    paste("CHHclass=CHHclass", levels(cF9$CHHclass), sep = "_"),
    unlist(sapply(h3ids, function(x) paste(x, "=", x, "_", levels(cF9[,x]), sep = ""))),
    paste(gsub("[-/\\?]", ".", supFac[-1]), "=", gsub("[-/\\?]", ".", supFac[-1]), "_TRUE", sep = "")
)

#rownames(pvalmatrix) <- sapply(strsplit(rownames(pvalmatrixsub),"="),function(x) x[2])

pvalmatrixsel <- pvalmatrix[selectP,]
pvalmatrixsel <- pvalmatrixsel[rowSums(abs(pvalmatrixsel) > 5) > 0,]

pvalmatrixsel <- pvalmatrixsel[-grep("(^tRNA=)|(^rRNA=)|(^pseudogene=)|(^miRNA=)|(^ncRNA=)|(^lincRNA=)", rownames(pvalmatrixsel)),]
pvalmatrixsub <- pvalmatrix[which(rowSums((pvalmatrix) > 5) > 1),]

cleanNames <- function(nammat) {    
    if(length(grep("FALSE", nammat)) > 0) pmat <- pmat[-grep("FALSE", nammat),]
    if(length(grep("overlaptype=", nammat)) > 0) pmat <- pmat[-grep("overlaptype=", nammat),]
    if(length(grep("not_known", nammat)) > 0) pmat <- pmat[-grep("not_known", nammat),]

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
    nammat <- gsub("med", "moderate", nammat)
    nammat <- gsub("isIR", "IR", nammat)
    nammat
}

cleanF9 <- cF9; colnames(cleanF9) <- cleanNames(colnames(cF9))
cleanF9 <- cbind(as.data.frame(gr9)[,1:3], cleanF9)
cleanF9 <- cleanF9[,-grep("mobile|mobType|AGO", colnames(cleanF9))]
cleanF9 <- cbind(LC = gr9$cluster, cleanF9)

write.table(cleanF9, col.names = NA, quote = FALSE, sep = "\t", file = "Loci_annotation.txt")

rownames(pvalmatrixselC) <- cleanNames(rownames(pvalmatrixsel))
rownames(pvalmatrixsubC) <- cleanNames(rownames(pvalmatrixsub))

pvalmatrixsubC <- pvalmatrixsubC[-grep("Ago", rownames(pvalmatrixsubC)),]
pvalmatrixsubC <- pvalmatrixsubC[-grep("mobType", rownames(pvalmatrixsubC)),]
pvalmatrixsubC <- pvalmatrixsubC[-grep("mobile", rownames(pvalmatrixsubC)),]

source("heatmap_centred.R")
pdf("~/My_Papers/loci_paper/figures/featureMatrix4_new_gr9.pdf",height = 12, width = 10)
heatmap.2(pvalmatrixsubC,
          Rowv = TRUE, Colv = NA, key = FALSE,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(6, 15), cexRow = 1.1, lhei = c(0.01, 5), lwid = c(0.5, 1))
dev.off()

pvalmatrixselC <- pvalmatrixselC[-grep("Ago", rownames(pvalmatrixselC)),]
pvalmatrixselC <- pvalmatrixselC[-grep("mobType", rownames(pvalmatrixselC)),]
pvalmatrixselC <- pvalmatrixselC[-grep("mobile", rownames(pvalmatrixselC)),]

pdf("~/My_Papers/loci_paper/figures/featureMatrix4b_gr9.pdf",width=8, height = 12)
heatmap.2(pvalmatrixselC,
          Colv = NA, Rowv = NA, key = FALSE,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(6, 15), cexRow = 1.1, lhei = c(0.01, 5), lwid = c(0.1, 1))
dev.off()

categories <- resMCA$desc.var$category
filcat <- lapply(categories, function(x) {
    x <- x[x[,5] > 0,]
    x <- x[x[,4] < 1e-4,]
    if(length(grep("overlaptype=", rownames(x))) > 0) x <- x[-grep("overlaptype=", rownames(x)),]
    if(length(grep("=.*NA", rownames(x))) > 0) x <- x[-grep("=.*NA", rownames(x)),]
    if(length(grep("=.*FALSE", rownames(x))) > 0) x <- x[-grep("=.*FALSE", rownames(x)),]
    if(length(grep("=.*not_known", rownames(x))) > 0) x <- x[-grep("=.*not_known", rownames(x)),]
    if(length(grep("=.*none", rownames(x))) > 0) x <- x[-grep("=.*none", rownames(x)),]
    if(length(grep("mobType", rownames(x))) > 0) x <- x[-grep("mobType", rownames(x)),]
    x
})

summariseType <- function(string, subset = 1:9) {
    z <- sapply(categories, function(x) x[rownames(x) == string,])
    z <- cbind(z, subsum = NA)
    z[1, 'subsum'] <- sum(z[1,subset], na.rm = TRUE)
    t(signif(z, 3))
}

summariseType("gene=gene_TRUE", 1:3)
summariseType("CHHclass=CHHclass_low", 1:3)
summariseType("tissue_expressionClass=specific", 1:3)
summariseType("easiRNA=easiRNA_TRUE", 4:9)
summariseType("CGclass=CGclass_high", 4:9)
summariseType("CHGclass=CHGclass_high", 4:9)
summariseType("CHHclass=CHHclass_high", 4:9)
summariseType("sizeclass=(0,50]", 1:3)
summariseType("promotor=promotor_TRUE")
summariseType("gene=gene_TRUE")
summariseType("hasH3K9Ac=hasH3K9Ac_low")
summariseType("NAT=NAT_TRUE")
summariseType("ratio_strand_class=strong bias")
summariseType("ratio21vs24Class=ratio21vs24Class_high")
summariseType("repetitivenessClass=repetitivenessClass_low")
summariseType("repetitivenessClass=repetitivenessClass_high")
summariseType("repetitivenessClass=repetitivenessClass_med")
summariseType("polIV_new=polIV_new_PolIV_independent")
summariseType("polV_new=polV_new_PolV_dependent", c(4,5,7))
summariseType("polV_new=polV_new_PolV_independent", c(6,8,9))
summariseType("hasH3K27me1=hasH3K27me1_high", c(4,5,7))
summariseType("tissue_expressionClass=inbetween", c(4,5,7))
summariseType("tissue_expressionClass=common", c(7:9))
summariseType("rdr2_new=rdr2_new_Rdr2_independent")
summariseType("predominant_5prime_letter=CGT", 8:9)
summariseType("ratio21vs24Class=ratio21vs24Class_low", c(4,5,7))
summariseType("RathE=RathE_TRUE", c(4,5,7))
summariseType("DNA=DNA_TRUE", c(6,8))
summariseType("mobType_mobileEffect.mobile.D3=mobType_mobileEffect.mobile.D3_TRUE", c(4,5,7))
summariseType("isIR=isIR_TRUE", c(6,8,9))
summariseType("hasH3K27me3=hasH3K27me3_med", c(6,8,9))
summariseType("hasH3K27me3=hasH3K27me3_low", c(6,8,9))
summariseType("hasH3K27me3=hasH3K27me3_high", c(5,7))
summariseType("hasH3K9me2=hasH3K9me2_med", c(6,8,9))
summariseType("hasH3K9me2=hasH3K9me2_low", c(6,8,9))
summariseType("hasH3K9me2=hasH3K9me2_low", c(6,8,9))
summariseType("hasH3K4me2=hasH3K4me2_med", 6:9)
summariseType("hasH3K4me3=hasH3K4me3_med", 6:9)
summariseType("hasH3K36me2=hasH3K36me2_med", 6:9)
summariseType("h3Class=h3Class_med", c(6,8,9))
summariseType("h3Class=h3Class_high", c(6,8,9))
summariseType("isPhased=isPhased_high", 7:9)
summariseType("isPhased=isPhased_moderate")
summariseType("isPhased=isPhased_none")
summariseType("dclTripleClass=dclTripleClass_TRUE", 7:9)
summariseType("CHG.mobileB=TRUE")
summariseType("CHG.mobileB=TRUE")
summariseType("hasExpression=hasExpression_TRUE")
summariseType("sizeclass=(2e+03,Inf]", 1:3)
summariseType("RC.Helitron=RC.Helitron_TRUE", c(4,5,7))
summariseType("LTR.Gypsy=LTR.Gypsy_TRUE", c(6,8,9))


lapply(1:nclust, function(ii) {
    x <- resMCA$desc.ind$para[[ii]]
    write.table(as.data.frame(gr9[as.integer(names(x)),])[,1:4], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA, file = paste("paragons_LC", ii, ".txt", sep = ""))
})

paralist <- list()
for(ii in 1:9) {
    paralist[[ii]] <- gr9[as.integer(names(resMCA$desc.ind$dist[[ii]]))]
    values(paralist[[ii]]) <- NULL
    values(paralist[[ii]])$cluster <- ii
}

paragons <- as.data.frame(do.call("c", paralist))
write.table(paragons, col.names = NA, file = "gr9_paragons.txt", quote = FALSE, sep = "\t")
