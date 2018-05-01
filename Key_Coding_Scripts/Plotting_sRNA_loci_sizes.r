#sRNA size graphs production
# /scripts/conscriptoR /home/bioinf/nem34/Plotting_sRNA_sizes.r
#Before removal of repetative data
#load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/aD_first_chlamy_segmentation_nick.RData")
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


#After removal of reptative data
load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/aDlt20_first_chlamy_segmentation_nick.RData")
#Before removal of repetative data
#load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/aD_first_chlamy_segmentation_nick.RData")

combined<-rbind(table(aD@alignments@ranges@width),table(rep(aD@alignments@ranges@width,aD@alignments$multireads)))
pdf("Figure2a.pdf",11,7)
#par(mfrow = c(2,1))
barplot(combined, beside = TRUE,legend.text=c("Non-Redundant","Redundant"),col=c("red","blue"),xlab="sRNA Size",ylab="Counts",main="sRNA Size Distribution - Non-Repetative")
dev.off()

pdf("Figure2b.pdf",11,7)
load("/home/bioinf/nem34/segmentMap_II/gr7_clustered.RData")
plot(density(log10(width(gr7))), main="Loci Size Distribution" , ylab="Density", xlab="Loci Size (log10)")
abline(v=log10(c(30,75,150,1000)))
dev.off()

#Before removal of repetative data

load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/aD_first_chlamy_segmentation_nick.RData")
combined<-rbind(table(aD@alignments@ranges@width),table(rep(aD@alignments@ranges@width,aD@alignments$multireads)))
pdf("sRNA_sizes_unfiltered.pdf")
barplot(combined, beside = TRUE,legend.text=c("Non-Redundant","Redundant"),col=c("red","blue"),xlab="sRNA Size",ylab="Counts",main="sRNA Size Distribution")
dev.off()






#gg <- ggplot(combined) + geom_density(aes(x=start),adjust=1/20, fill="red") + facet_grid(annot~chrom, scales = "free") +  theme_bw()
#ggsave(gg,file=paste0(prefix,"Clustercoverage.eps"),width=30,height=5)
#ggsave(gg,file=paste0(prefix,"Clustercoverage.png"),width=30,height=5)

#non-redundant
#pdf()
#plot(table(aD@alignments@ranges@width))
#dev.off()
#height<-rbind(table(aD@alignments@ranges@width), table(rep(aD@alignments@ranges@width,aD@alignments$multireads)))
#pdf("sRNA_sizes_lt20.pdf")
# barplot(height, beside = TRUE)
#dev.off()
#table(rep(aD@alignments@ranges@width,aD@alignments$multireads))


#pdf("sRNA_sizes_lt20.pdf")
#barplot(table(aD@alignments@ranges@width), col="red")
#par(new=T)
#barplot(table(rep(aD@alignments@ranges@width,aD@alignments$multireads)), col="green")
#dev.off()
