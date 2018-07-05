library(baySeq)
library(segmentSeq)

cl <- makeCluster(24)

setwd("~/Code/segmentMap_III")

chrlens <- c(30427671, 19698289, 23459830, 18585056, 26975502)

samples <- read.delim("Arabidopsis_samples.txt", as.is = TRUE)

samples$files[samples$Library == "SL12"] <- "/data/pipeline/prod/SL12/SL12.ID15_FC5372.lane3.reads.filtered_7_5_2008.v_Arabidopsis_thaliana_genome-tair9.patman.gff2"
samples$files[samples$Library == "SL9"] <- "/data/pipeline/prod/SL9/SL9.ID15_FC5372.lane1.reads.filtered_7_5_2008.v_Arabidopsis_thaliana_genome-tair9.patman.gff2"
samples$files[samples$Library == "SL11"] <- "/data/pipeline/prod/SL11/SL11.ID15_FC5372.lane2.reads.filtered_7_5_2008.v_Arabidopsis_thaliana_genome-tair9.patman.gff2"
samples$files[samples$Library == "SL4_1"] <- "/data/pipeline/prod/SL4/SL4.FC5365.combined.reads.filtered_7_5_2008.v_Arabidopsis_thaliana_genome-tair9.patman.gff2"

mirFiles <- gsub("genome-tair9.patman.gff2", "mirbase-hairpin.patman.gff2", samples$files)
mirFiles <- mirFiles[grep("mirbase-hairpin.patman.gff2", mirFiles)]
mirTags <- unique(unlist(sapply(mirFiles, function(x) {
  tryTag <- try(unique(read.table(x, header = FALSE, as.is = TRUE)[,10]))
  if(class(tryTag) == "try-error") return(NULL) else return(tryTag)
})))

trnaFiles <- gsub("genome-tair9.patman.gff2", "trna.patman.gff2", files)
trnaFiles <- trnaFiles[grep("trna.patman.gff2", trnaFiles)]
trnaTags <- unique(unlist(lapply(trnaFiles, function(x) {
  tryTag <- try(unique(read.table(x, header = FALSE, as.is = TRUE)[,10]))
  if(class(tryTag) == "try-error") return(NULL) else return(tryTag)
})))

rrnaFiles <- gsub("genome-tair9.patman.gff2", "rrna.patman.gff2", files)
rrnaFiles <- rrnaFiles[grep("rrna.patman.gff2", rrnaFiles)]
rrnaTags <- unique(unlist(lapply(rrnaFiles, function(x) {
  tryTag <- try(unique(read.table(x, header = FALSE, as.is = TRUE)[,10]))
  if(class(tryTag) == "try-error") return(NULL) else return(tryTag)
})))

rawFiles <- paste("/data/pipeline/prod/", samples$Library, "/", samples$Library, ".filtered_trimmed_reads.fasta", sep = "")
rawFiles[samples$Owner == "Sebastian"] <- dir("/home/sm934/seqdata/sRNAexternallibs", full.name = TRUE, pattern = "nonredundant\\.fa$")

samples$Library[is.na(samples$Library)] <- gsub("_nonredundant.aligned.*", "", gsub(".*/", "", samples$files[is.na(samples$Library)]))

for(ii in 1:length(rawFiles))
  if(!file.exists(paste("~/Code/segmentMap_II/alignJunk/", samples$Library[ii], "_scored_patman.out", sep = "")))
    system(paste("perl /scripts/run_standard_patman.pl -D /home/tjh48/Code/arabidopsisFilter/arabidopsisTAIR9_filterJunk.fasta -P ", rawFiles[ii], " -o ~/Code/segmentMap_II/alignJunk -n ", samples$Library[ii], sep = ""))

system("perl /scripts/run_standard_patman.pl -D /home/tjh48/Code/arabidopsisFilter/arabidopsisTAIR9_filterJunk.fasta -P /data/pipeline/prod/SL9/SL9.ID15_FC5372.lane1.reads.filtered_7_5_2008.fasta -o alignJunk -n SL9")
system("perl /scripts/run_standard_patman.pl -D /home/tjh48/Code/arabidopsisFilter/arabidopsisTAIR9_filterJunk.fasta -P /data/pipeline/prod/SL11/SL11.ID15_FC5372.lane2.reads.filtered_7_5_2008.fasta -o alignJunk -n SL11")
system("perl /scripts/run_standard_patman.pl -D /home/tjh48/Code/arabidopsisFilter/arabidopsisTAIR9_filterJunk.fasta -P /data/pipeline/prod/SL12/SL12.ID15_FC5372.lane3.reads.filtered_7_5_2008.fasta -o alignJunk -n SL12")


junkTags <- unique(unlist(lapply(dir("~/Code/segmentMap_II/alignJunk", pattern = "scored_patman.out", full.names = TRUE), function(x) unique(read.delim(x, header = FALSE, as.is = TRUE)[,2]))))

failed <- c("SL27", "SL33")
garbage <- c("SL44_2", "SL103_1", "SL53")

files <- extant$files[1:84]
dirs <- gsub("(.*/)SL.*", "\\1", files)
bamfile <- sapply(dirs, function(dir) dir(dir, pattern = "genome-tair9.patman.bam$", full.names = TRUE))

for(bf in bamfile)
    system(paste('lftp -u tjh48@cam.ac.uk,fire5fly -e "mput ', bf, ';quit" webin.ebi.ac.uk', sep = ""))
sf <- samples[samples$files %in% files,]
sfdf <- data.frame(sample_alias = sf[,3], tax_id = 3702, scientific_name = "Arabidopsis thaliana", common_name = "thale cress",
           sample_title = sf[,3], sample_description = sf[,2], tissue_type = gsub(" .*", "", sf[,1]), replicate_group = sf$Replicate.group, technical_group = sf$Technical.Group)              
seqCentre <- read.csv("seq_centre.csv")

sfdf <- cbind(sfdf, seqCentre[match(sfdf[,1], gsub("RUN_", "", seqCentre[,1])),2])
sfdf <- cbind(sfdf, bamfile, md5)
sfdf[,11] <- gsub(".*/", "", sfdf[,11])
write.table(sfdf, file = "~/sfdf.txt", sep = "\t", quote = FALSE, row.names = FALSE)


extant <- samples[!(samples$Library %in% c(failed, garbage)),]
readFiles <- samples$files[!(samples$Library %in% c(failed, garbage))]

#readFiles <- readFiles[1:4]
#extant <- extant[1:4,]

aD <- readGeneric(files = readFiles, dir = "", replicates = extant$Replicate.group, libnames = extant$Library, chrs = c("1", "2", "3", "4", "5"), gap = 100, chrlens = chrlens, header = FALSE, cols = c(chr = 1, start = 4, end = 5, count = 6, strand = 7, tag = 10), polyLength = 10)

save(aD, file = "aD.RData")
load("aD.RData")

aD <- aD[!(as.character(values(aD@alignments)$tag) %in% junkTags),]

save(aD, file = "aD_noJunk.RData")
load("aD_noJunk.RData")

unqDat <- aD[!duplicated(aD@alignments$tag),]

lenDat <- lapply(15:29, function(bp)
    colSums(unqDat@data[width(unqDat@alignments) == bp,]))
lenDat <- do.call("cbind", lenDat)
colnames(lenDat) <- 15:29

lenDat / rowSums(lenDat)
plot(NA, NA, xlim = c(15,29), ylim = c(0,1), xlab = "sRNA length (bp)", ylab = "Proportion of reads")
redblue <- colorRampPalette(c("red", "blue"))(100)
apply(lenDat / rowSums(lenDat), 1, function(y) lines(y, x = 15:29, col = redblue[round(100 * y["21"] / (y["21"] + y["24"]))]))
savePlot("length_distributions.png")


sD <- processAD(aD, gap = NULL, verbose = TRUE, cl = cl)

save(sD, file = "sD.RData")
load("sD.RData")

cutoffSeg <- heuristicSeg(sD = sD, RKPM = 1000, gap = 100, getLikes = FALSE, cl = cl)
save(cutoffSeg, file = "cutoffSeg.RData")

load("aD_noJunk.RData")
load("cutoffSeg.RData")

gc()

cD <- cutoffSeg; newCounts = FALSE; bootStraps = 3; inferNulls = TRUE; nasZero = FALSE; usePosteriors = TRUE

cutoffSegLike <- lociLikelihoods(cutoffSeg, aD, newCounts = FALSE, bootStraps = 3, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl = cl)

save(cutoffSegLike, file = "cutoffSegLike.RData")
rm(aD)
gc()

load("sD.RData")

load("cutoffSegLike.RData")
cD <- cutoffSegLike
lociCutoff = 0.5; nullCutoff = 0.5; subRegion = NULL; getLikes = FALSE; lR = FALSE; samplesize = 1e5; cl = cl; tempDir = "tmpFile"; subRegion = NULL#subRegion = data.frame(chr = "1", start = 1, end = 1e6)
largeness <- 1e8

classSeg <- classifySeg(sD = sD, cD = cD, aD = aD, 
                        lociCutoff = 0.9, nullCutoff = 0.9, getLikes = FALSE, lR = FALSE, samplesize = 1e5, cl = cl, tempDir = "tmpFile2")

pdf("~/deleteMe.pdf")
plotGenome(aD = aD, classSeg, chr = "1", limits = c(0, 10000))
dev.off()

save(classSeg, file = "new_classSeg.RData")

load(file = "new_classSeg.RData")

classSegLike <- lociLikelihoods(
                  cD = classSeg, aD, newCounts = FALSE, bootStraps = 3, inferNulls = TRUE, nasZero = FALSE, usePosteriors = TRUE, cl = cl
                  )
save(classSegLike, file = "new_classSegLike.RData")

load("new_classSegLike.RData")

png("exampleSeg.png", height = 800, width = 600)
plotGenome(aD = aD, classSegLike, chr = "1", limits = c(16957, 17154) + c(-10000, 10000), showNumber = FALSE)
dev.off()


load("classSegLike.RData")
genome(classSegLike@coordinates) <- NA

splitNulls <- lapply(1:length(dir("tmpFile", pattern = "subNull")), function(ii) {
  load(file = paste(tempDir, "/subNull_", ii, ".RData", sep = ""))
  emptyD <- rowSums(sapply(1:ncol(potnullD), function(jj) as.integer(potnullD@data[,jj]))) == 0
  values(potnullD@coordinates)$empty <- FALSE
  values(potnullD@coordinates)$empty[emptyD] <- TRUE
  potnullD@data <- DataFrame(matrix(ncol = length(potnullD@replicates), nrow = 0))
  potnullD <- potnullD[rowSums(sapply(1:ncol(potnullD@locLikelihoods), function(ii) potnullD@locLikelihoods[,ii] > -Inf), na.rm = TRUE) > 0 | emptyD,]
  potnullD
})

save(splitNulls, file = paste(tempDir, "/splitNulls.RData", sep = ""))

load("tmpFile/splitNulls.RData")

sD@locLikelihoods <- DataFrame(do.call("rbind", lapply(1:length(dir(tempDir, pattern = "subLoc")), function(ii) {
  load(file = paste(tempDir, "/subLoc_", ii, ".RData", sep = ""))
  log(subLoc)
})))

save(sD, file = "tmpFile/sD.RData")

nullD <- .mergeSegData(splitNulls)

save(nullD, file = "tmpFile/nullD.RData")

locNums <- cbind(replicates = levels(classSegLike@replicates), extant[match(levels(classSegLike@replicates), extant$Replicate.group),c(1,2,4)], numberOfLoci = colSums(exp(classSegLike@locLikelihoods) > 0.5))
locNums <- locNums[order(locNums$numberOfLoci),]

locNums$Libraries <- sapply(locNums$replicates, function(rep) paste(extant$Library[extant$Replicate.group == rep], sep = ",", collapse = ","))
locNums$libsizes <- sapply(locNums$replicates, function(rep) mean(classSegLike@libsizes[classSegLike@replicates == rep]))

corLoc <- matrix(NA, nrow = ncol(classSegLike@locLikelihoods), ncol = ncol(classSegLike@locLikelihoods))
for(ii in 1:ncol(classSegLike@locLikelihoods))
  for(jj in 1:ncol(classSegLike@locLikelihoods))
  corLoc[ii,jj] <- cor(exp(classSegLike@locLikelihoods[,ii]), exp(classSegLike@locLikelihoods[,jj]))
rownames(corLoc)[locNums$replicates] <- colnames(corLoc)[locNums$replicates] <- paste(locNums$Tissue.type, locNums$Plant.Expt.Type, locNums$Owner, locNums$Libraries, sep = "_")

library(dendrogram)
locLibsizes <- locNums$libsizes[match(1:nrow(corLoc), locNums[,1])]
numLocs <- locNums$numberOfLoci[match(1:nrow(corLoc), locNums[,1])]

plot(locLibsizes, colSums((classSegLike@locLikelihoods > log(0.5))[rowSums(classSegLike@locLikelihoods > log(0.5)) == 1,]))

exp(classSegLike@locLikelihoods[,1])

par(mar = c(5,4,4,15))
rc <- rgb((locLibsizes / max(locLibsizes))^0.5, (locLibsizes / max(locLibsizes))^0.5, (locLibsizes / max(locLibsizes))^0.5)
cc <- rc <- rgb((numLocs / max(numLocs))^0.5, (numLocs / max(numLocs))^0.5, (numLocs / max(numLocs))^0.5)
heatmap(corLoc, labCol = NA, RowSideColors = rc, ColSideColors = cc, symm = TRUE)

png("locNum_libsizes.png", height = 600, width = 600)
plot(y = locNums$numberOfLoci, x = locNums$libsizes, xlab = "Library sizes", ylab = "Number of loci")
dev.off()

annotation <- read.delim("/data/public_data/arabidopsis/TAIR_10/TAIR10_GFF3_genes_transposons.gff", header = FALSE, as.is = TRUE)
annotation[annotation[,7] == ".",7] <- "*"
annotationGR <- GRanges(seqnames = gsub("Chr", "", annotation[,1]), IRanges(start = annotation[,4], end = annotation[,5]), strand = annotation[,7])

segOver <- getOverlaps(classSegLike@coordinates, annotationGR, whichOverlaps = TRUE, cl = NULL)

locOver <- segOver[which(rowSums(exp(classSegLike@locLikelihoods) > 0.5) > 0)]
loci <- classSegLike[rowSums(classSegLike@locLikelihoods > log(0.5)) > 0,]

checkAnn <- function(overlaps, annType)
  any(annType %in% annotation[overlaps,3])

locann <- sapply(unique(annotation[,3]), function(annType) {
  sapply(locOver, checkAnn, annType = annType)
})

annNumbers <- do.call("cbind", lapply(1:ncol(loci@locLikelihoods), function(ii)
                                      colSums(locann[loci@locLikelihoods[,ii] > log(0.5),])
                                      ))

unqNumbers <- do.call("cbind", lapply(1:ncol(loci@locLikelihoods), function(ii)
                                      colSums(locann[loci@locLikelihoods[,ii] > log(0.5) & rowSums(loci@locLikelihoods > log(0.5)) == 1,,drop = FALSE])
                                      ))

repNames <- extant[match(levels(loci@replicates), extant[,6]),]
repNames <- paste(repNames$Tissue.type, repNames$Plant.Expt.Type, repNames$Owner, repNames$Libraries, sep = "_")

colnames(unqNumbers) <- colnames(annNumbers) <- repNames
normAnn <- t(t(annNumbers) / annNumbers[1,])[-1,]

png("annotation.png", height = 600, width = 1000)
par(mar = c(1, 15, 1, 1))
barplot(normAnn[,as.numeric(as.character(locNums[order(locNums$Owner),1]))], horiz = TRUE, las = 2, col = rainbow(20), legend =TRUE, cex.names = 0.8)
barplot(normAnn[,order(normAnn[1,], decreasing = TRUE)], horiz = TRUE, las = 2, col = rainbow(20), legend =TRUE, cex.names = 0.8)
barplot(unqNumbers[,order(unqNumbers[1,], decreasing = TRUE)], horiz = TRUE, las = 2, col = rainbow(20), legend =TRUE, cex.names = 0.8)
dev.off()

sapply(1:length(levels(loci@replicates)), function(ii) {
  write.table(as.data.frame(loci@coordinates[loci@locLikelihoods[,ii] > log(0.5),]), file = paste("loci_coordinates/", gsub(" ", "", gsub("/", "_plus_", gsub("_$", ".txt", repNames[ii]))), sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
  })

sapply(unique(annotation[,3]), function(annType) {
  redType <- annotation[annotation[,3] == annType,c(1, 4, 5)]
  colnames(redType) <- c("seqnames", "start", "end")
  redType[,1] <- gsub("Chr", "", redType[,1])
  write.table(redType, file = paste("ann_coordinates/", annType, ".txt", sep = ""), sep = "\t", quote = FALSE, row.names = FALSE)
})

bb <- lapply(unique(annotation[,3]), function(annType) {
  bbSamp <- lapply(1:length(repNames), function(ii) {
    blockBootstrap(p
                   file1 = paste("loci_coordinates/", gsub(" ", "", gsub("/", "_plus_", gsub("_$", ".txt", repNames[ii]))), sep = ""),
                   file2 = paste("ann_coordinates/", annType, ".txt", sep = ""),
                   chrLengths = chrLengths)
  })
  names(bbSamp) <- repNames
  bbSamp
})

names(bb) <- unique(annotation[,3])

x11()
pdf("Sebastian_funny_loci.pdf", height = 20, width = 30)
plotGenome(aD, classSegLike, chr = "1", limits = c(4987695, 5009703) + c(-500, 500))
plotGenome(aD, classSegLike, chr = "1", limits = c(4987695, 5009703) + c(-500, 500))
plotGenome(aD, classSegLike, chr = "1", limits = c(5356765,5357154) + c(-500, 500))
dev.off()
