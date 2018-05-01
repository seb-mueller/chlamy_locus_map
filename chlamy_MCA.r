#Adapted by Nick Matthews from Tom Hardcastle's orginal code
#Date: 01/03/16

#To submit this script to condor
#		 /scripts/conscriptoR /home/bioinf/nem34/segmentMap_II/chlamy_MCA_8.r -p19

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

source("/home/sm934/code/R/seb_functions.r")
load("loci7_fdr01_c.rdata") #not mentioned, should it be loci7?
load("gr7_all_annotations15.rdata")

##MCA for categoriacal data!
#multivariate methods that allows us to analyze the systematic patterns of variations with categorical data
#keep in mind that MCA applies to tables in which the observations are described by a set of qualitative (i.e. categorical) variables

# selected factors which will be used to inform the clustering
selFac <- c("sizeclass","predominant_5prime_letter","ratio21vs20Class","expressionClass", 
	"repetitivenessClass","Phased","meth","zygotespecific","vegetativespecific",
	"CC125specific","CC1883specific","notDCL3dependent","DCL3dependent",
	"miRNA","gene","IR","TR","promotor","ratio_strand_class","SoloLTR","TEorder") 

# supplementary factors for which association with clusters will be calculated, but which will not inform the clustering
supFac <- c("CDS","exon","fiveprimeUTR","threeprimeUTR","Jspecific","TE","TEclass","methCG","methCHH","methCHG","ratioBigvsNormalClass","ratioSmallvsNormalClass")

cF7 <- as.data.frame(elementMetadata(gr7[,c(selFac,supFac)]))                                                      

png("catchall.png")

#tables summarising output
write("", file = "ClassTable_gr7.txt")
for(ii in 1:ncol(cF7)) {
    cat(colnames(cF7)[ii], "\t", paste(levels(cF7[,ii]), collapse = "\t"), "\n", file = "ClassTable_gr7.txt", append = TRUE)
    cat("", "\t", paste(as.numeric(table(cF7[,ii])), collapse = "\t"), "\n\n", file = "ClassTable_gr7.txt", append = TRUE)
}

# MCA
mc7 <- MCA(cF7, graph = FALSE,ncp = 4 ,quali.sup=which(colnames(cF7) %in% supFac))

# parameter sweep on dimensions 1-15 and clusters 2-15
cl <- makeCluster(14)
dimList <- list()
for(nn in 1:15) {
    mc7 <- MCA(cF7, graph = FALSE,ncp = nn ,quali.sup=which(colnames(cF7) %in% supFac))
    dimList[[nn]] <- c(list(NA), parLapply(cl, 2:15, function(i, coords)
        kmeans(x = coords, iter.max = 1000, nstart = 1000, centers = i), coords = mc7$ind$coord))
}
    
save(dimList, file = "dimList.RData")
#load("dimList.RData")

# probably not useful as a measure of performance
png("kmean_ss.png")
image(sapply(dimList, function(x) sapply(x[-1], function(y) y$betweenss / y$totss)))
dev.off()

# rand test to compare overlap with transposable element superfamilies 
zzz <- (sapply(dimList, function(x) sapply(x[-1], function(y) adjustedRandIndex(y$cluster, gr7$TE))))
png("kmean_randTE.png")
image(zzz)
dev.off()

# rand test to compare overlap with annotated features. Used in Tom's code, I no longer use Overlaptype
#zz <- (sapply(dimList, function(x) sapply(x[-1], function(y) adjustedRandIndex(y$cluster, gr7$overlaptype)))) #Overlaptype loses lots of important info
#png("kmean_rand.png")
#image(zz)
#dev.off()

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

save(dimStab, file = "dimStab.RData")
#load("dimStab.RData")

#Plot stability analysis

pdf("stabilityplots.pdf",20,20)
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

pdf("davies_bouldin_image.pdf")
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

pdf("dunn_image.pdf")
image(dunnM)
image(dunnM2)
image(dunnM3)
dev.off()

#Examine images and choose correct dimension/cluster number. Need to talk to Tom for interpretation


# HCPC code from FactoMiner needs tweak to work on kmeans only.
source("/home/tjh48/Code/segmentMap_III/hcpc.R")

#MCA with clusters and dimensions set according to images
nclust <- 4; ndim <- 3
#nclust <- 4; ndim <- 4
#nclust <- 5; ndim <- 4
#nclust <- 6; ndim <- 4
mc7 <- MCA(cF7, graph = FALSE,ncp = ndim ,quali.sup=which(colnames(cF7) %in% supFac))
resMCA <- HCPC(mc7, graph = FALSE, proba = 1, consol = FALSE, order = FALSE, nb.clust = nclust, kk = nclust, method = "centroid")
gr7$cluster <- as.factor(resMCA$data.clust$clust)

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


colnames(pvalmatrix) <- paste(colnames(pvalmatrix), " (", sapply(colnames(pvalmatrix), function(ii) sum(gr7$cluster == ii)), ")", sep = "")

# select only significant classes
pvalmatrixsub <- pvalmatrix[which(rowSums(abs(pvalmatrix) > 3) > 1),]
pvalmatrixsub <- pvalmatrixsub[-grep("FALSE", rownames(pvalmatrixsub)),]
#pvalmatrixsub <- pvalmatrixsub[-grep("overlaptype=", rownames(pvalmatrixsub)),]

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
	"CDS=CDS_TRUE","exon=exon_TRUE","fiveprimeUTR=fiveprimeUTR_TRUE","threeprimeUTR=threeprimeUTR_TRUE","Jspecific=Jspecific_TRUE",
	"miRNA=miRNA_TRUE",
	"IR=IR_TRUE",
	"TR=TR_TRUE",
	paste("TE=TE_", c("none","LTR.Copia","LTR.Gypsy","LTR.TOC1","LTR.REM1","DIRS","LINE.RTE","LINE.L1","LINE.Dualen/RandI","SINE",
						"TIR.TCR1","TIR.TOC2","TIR.hAT","TIR.Novosib","TIR.Gulliver","TIR.Mariner","TIR.P","TIR.EnSpm",
						"DNA.NonAut","Un.TE1","Un.TE2","NonLTR") , sep = ""),
	paste("TEorder=TEorder_", c("LTR","DIRS","LINE","SINE","TIR") , sep = ""),
	paste("TEclass=TEclass_", c("DNA","RET") , sep = "")
	)
	

pvalmatrixsel <- pvalmatrix[selectP,]

# heatmap function needs some tweaks to work on these data
source("/home/bioinf/tjh48/Code/segmentMap_III/heatmap_centred.R")


#examine heatmaps
for(ii in 1:nrow(pvalmatrixsub))
    rownames(pvalmatrixsub)[ii] <- gsub(paste("=", gsub("=.*", "", rownames(pvalmatrixsub)[ii]), "_", sep = ""), "=", rownames(pvalmatrixsub)[ii])
png("featureMatrix4_new_gr7.png",width=1200, height = 1600)
par(mar = c(5.1, 4.1, 9.1, 2.1))
heatmap.2(pvalmatrixsub,
          Colv = NA,
          col = rev(c(colorRampPalette(colors = c("blue", "white"))(255), colorRampPalette(colors = c("white", "red"))(255))), scale = "none", margins = c(8, 20), cexRow = 1.1)
dev.off()

for(ii in 1:nrow(pvalmatrixsel))
    rownames(pvalmatrixsel)[ii] <- gsub(paste("=", gsub("=.*", "", rownames(pvalmatrixsel)[ii]), "_", sep = ""), "=", rownames(pvalmatrixsel)[ii])
png("featureMatrix4b_gr7.png",width=1200, height = 1600)
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

save(gr7, file="gr7_clustered.RData")
save(resMCA, file="resMCA.RData")
#load("resMCA.RData")


# try plotting density of loci in different clusters across genome; compare with gene and transposable element densities.

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

#New methylation
methCG<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CGmeth.gff3")
annottrack_methCG <- data.frame(chrom=as.factor(as.character(methCG@seqnames)), start=start(methCG), annot=rep("meth2",length(methCG)))
methCHH<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHHmeth.gff3")
annottrack_methCHH <- data.frame(chrom=as.factor(as.character(methCHH@seqnames)), start=start(methCHH), annot=rep("meth2",length(methCHH)))
methCHG<-import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHGmeth.gff3")
annottrack_methCHG <- data.frame(chrom=as.factor(as.character(methCHG@seqnames)), start=start(methCHG), annot=rep("meth2",length(methCHG)))
annottrack_meth2 <- rbind(annottrack_methCG,annottrack_methCHH,annottrack_methCHG)

#loci
gr<-gr7;prefix="loci_"
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



lapply(1:nclust, function(ii) {
    x <- resMCA$desc.ind$para[[ii]]
    write.table(as.data.frame(gr7[as.integer(names(x)),])[,1:4], sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA, file = paste("paragons_LC", ii, ".txt", sep = ""))
})

dev.off() 