histoneOverlap <- function(loci, histoneFile = "/projects/tjh48/histone/histoneGR.RData")
    {
        load(histoneFile)
        peakOverlaps <- do.call("cbind.data.frame", lapply(peakGR, function(peaks) {
            getOverlaps(loci, peaks, whichOverlaps = FALSE, ignoreStrand = FALSE, cl = NULL)
        }))
        values(loci) <- cbind(values(loci), peakOverlaps)
        loci
    }
        
    
selectLoci <- function (cD, likelihood, FDR, FWER, perReplicate = TRUE, returnBool = FALSE) 
{
    
    if (!missing(likelihood)) {
        selLoc <- cD@locLikelihoods > log(likelihood)
        if (returnBool) 
            return(selLoc)
        else selLoc <- which(rowSums(selLoc) > 0)
    }
    else {
        if (!missing(FDR)) {
            controlFunction <- segmentSeq:::.controlFDR
            controlCrit <- FDR
        }
        else if (!missing(FWER)) {
            controlFunction <- segmentSeq:::.controlFWER
            controlCrit <- FWER
        }
        else stop("No criterion for locus selection given.")
        if (perReplicate) {
            selRep <- lapply(1:ncol(cD@locLikelihoods), function(jj) controlFunction(exp(cD@locLikelihoods[, 
                jj]), controlCrit))
            if (returnBool) {
                bool <- do.call("cbind", lapply(1:length(selRep), 
                  function(ii) {
                    selBool <- rep(FALSE, nrow(cD))
                    if (length(selRep[[ii]]) > 0) 
                      selBool[selRep[[ii]]] <- TRUE
                    selBool
                  }))
                colnames(bool) <- colnames(cD@locLikelihoods)
                return(bool)
            }
            selLoc <- sort(unique(unlist(selRep)))
        }
        else {
            selLoc <- controlFunction(1 - exp(rowSums(log(1 - 
                exp(cD@locLikelihoods)))), controlCrit)
            if (returnBool) {
                bool <- rep(FALSE, nrow(cD))
                bool[selLoc] <- TRUE
                return(bool)
            }
        }
    }
    if (length(selLoc) == 0) 
        stop("No loci found for the given selection criterion.")
    cD[selLoc, ]
}

easiRNAS <- function(locAnn) {
    esd <- read.csv("nature13069-s1.csv", skip = 1)
    grESD <- GRanges(seqnames = gsub("Chr", "", esd[,1]), IRanges(start = esd[,2], end = esd[,3]))
    gr7$easiRNA <- getOverlaps(gr7, grESD, whichOverlaps = FALSE)
    gr7
}

mobileMethClasses <- function(locAnn) {
    load("/home/tjh48/Code/root_grafts_methylation/newFigures/srnaMethOverlaps_newLoci.RData")
    cgMobClass <- do.call("cbind", lapply(cgovCods, function(x) {
        seqlevels(x) <- gsub("Chr", "", seqlevels(x))
        getOverlaps(gr7, x, whichOverlaps = FALSE)
    }))
    colnames(cgMobClass) <- paste("CG.mobile", colnames(cgMobClass), sep = "")

    chgMobClass <- do.call("cbind", lapply(chgovCods, function(x) {
        seqlevels(x) <- gsub("Chr", "", seqlevels(x))
        getOverlaps(gr7, x, whichOverlaps = FALSE)
    }))
    colnames(chgMobClass) <- paste("CHG.mobile", colnames(chgMobClass), sep = "")
    
    chhMobClass <- do.call("cbind", lapply(chhovCods, function(x) {
        seqlevels(x) <- gsub("Chr", "", seqlevels(x))
        getOverlaps(gr7, x, whichOverlaps = FALSE)
    }))
    colnames(chhMobClass) <- paste("CHH.mobile", colnames(chhMobClass), sep = "")

    values(gr7) <- cbind(values(gr7), DataFrame(cgMobClass), DataFrame(chgMobClass), DataFrame(chhMobClass))

    gr7
}


phaseMatch <- function(locAnn)
    {
        ##TASI analyse

        beds <- as.data.frame(locAnn)[,1:3]
        beds[,1] <- paste("Chr", beds[,1], sep = "")
        rownames(beds) <- 1:nrow(beds)
        write.table(beds, col.names = TRUE, row.names = TRUE, file = "phaser/phasloc.txt", quote = FALSE, sep = "\t")
        system("python phaser/run_phasing.py")
        
        tasi <- read.csv("phaser/phasing_results_by_locus_21nt.tsv",sep="\t",header=FALSE)
                                        #table(locAnn$cluster[tasi[,"V10"]< c(-5)] & tasi[,"V4"]< c(-5)])
        isphased <- rep("none",nrow(tasi))
        tasi[,4:19] <- exp(tasi[,4:19])
        tasi[,4:19] <- apply(tasi[,4:19], 2, function(x) {x[x<1] <- p.adjust(x[x < 1], method = "BH"); x})
        isphased[rowSums(tasi[,4:19] < 0.05)>1] <- "moderate"
        isphased[rowSums(tasi[,4:19] < 0.05)>5] <- "high"
        
        tasi[,1] <- gsub("Chr", "", tasi[,1])
        tasmat <- match(apply(as.data.frame(locAnn)[,1:3], 1, function(x) paste(gsub(" ", "", x), collapse = ":")),
                        apply(as.data.frame(tasi)[,1:3], 1, function(x) paste(gsub(" ", "", x), collapse = ":")))
        locAnn$isPhased <- "none"
        locAnn$isPhased[!is.na(tasmat)] <- isphased[tasmat[!is.na(tasmat)]]        
        locAnn$isPhased <- ordered(as.factor(locAnn$isPhased), levels = c("none", "moderate", "high"))
        
        locAnn
    }

TIGR <- function()
    {
        transposonsgff <- import.gff3("TAIR10_Transposable_Elements.gff") #GRanges
        transposonsgff$type <- as.factor(transposonsgff$Super_Family); tmp <- transposonsgff$Name;transposonsgff$Name <- NULL
        transposonsgff$ID <- transposonsgff$Family;  transposonsgff$Name <- tmp
        transposonsgff$Family <- NULL;transposonsgff$Super_Family <- NULL
        
        transposonsgff$type <- combine_factor(transposonsgff$type, c(1:12,13,13,13,16:18))
        levels(transposonsgff$type)[13] <- "RathE"
        
                                        #lincRNAs
        lincRNAs <- import.gff3("/home/sm934/seqdata/lincRNAs_Liu/LincRNAs_identified_by_both_RepTAS_and_RNA-seq.gff")
        lincRNAs$Name <- lincRNAs$ID
        levels(lincRNAs$type) <- "lincRNA"
        
                                        #natRNAs
        natRNAs <- import.gff3("/home/sm934/seqdata/Luo2013_NATs_stringentb.gff3")[,1:6]
        natRNAs$type <- "NAT"
        
                                        #inverted repeats from irfinder
                                        #wine irf305.dos.exe arabidopsis.fa 2 3 5 80 10 40 500000 10000 -d -h -t4 74 -t5 493 -t7 10000
                                        #http://tandem.bu.edu/irf/irf.unix.help.html
        irs <- import.gff3("/data/public_data/arabidopsis/TAIR_10/custom_tracks/arabidopsis.fa.2.3.5.80.10.40.500000.10000_b.gff3")
        
                                        #/data/public_data/arabidopsis/TAIR_10/TAIR10_GFF3_genes_transposons.gff
        anno <- import.gff3("/data/public_data/arabidopsis/TAIR_10/TAIR10_GFF3_genes_transposons.gff")
                                        #annogr <- GRanges(seqnames=anno$space,ranges = IRanges(start=start(anno),end=end(anno)),strand=strand(anno),type=anno$type,group=anno$group)
        annogr <- anno[!anno@elementMetadata$type=="chromosome",]
                                        #transposons:
                                        #ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt
        
        miRNAs <- annogr[annogr@elementMetadata$type=="miRNA",]
                                        #transposons <- annogr[annogr@elementMetadata$type=="transposable_element" | annogr@elementMetadata$type=="transposable_element_gene" | annogr@elementMetadata$type=="transposon_fragment",]
                                        #there is a better annotation on TAIR10 ftp://ftp.arabidopsis.org/home/tair/Genes/TAIR10_genome_release/TAIR10_transposable_elements/TAIR10_Transposable_Elements.txt
        
        ncRNA <- annogr[annogr@elementMetadata$type=="ncRNA",]
        mRNA <- annogr[annogr@elementMetadata$type=="mRNA",]
        elements <- annogr[annogr@elementMetadata$type=="miRNA" |annogr@elementMetadata$type=="ncRNA" | annogr@elementMetadata$type=="gene" | annogr@elementMetadata$type=="pseudogene" | annogr@elementMetadata$type=="rRNA" |  annogr@elementMetadata$type=="tRNA",] #only 11 sno,sn RNA, therefore discarded
        elements$Index <- NULL;elements$Note <- NULL;elements$Parent <- NULL;elements$Index <- NULL;elements$Derives_from <- NULL;elements$Alias <- NULL
        rRNA <- annogr[annogr@elementMetadata$type=="rRNA",]
        tRNA <- annogr[annogr@elementMetadata$type=="tRNA",]
        
        
        promotors1k <- flank(elements[elements$type=="gene"],500)
        promotors1k$type <- droplevels(promotors1k$type)
        levels(promotors1k$type) <- "promotor"
        
        genes <- c(elements,transposonsgff,lincRNAs,promotors1k,natRNAs)
        genes$type <- droplevels(genes$type)
        genes$sizeclass <- cut(width(genes),breaks=c(0,500,2000,max(width(genes))))
        levels(genes$sizeclass) <- c("Short","Medium","Long")
        genes$isTE <- genes$type %in% levels(transposonsgff$type)
        
#        seqsmRNA <- DNAStringSet(vector(mode = "character", length = length(mRNA)))
#        100mb and therefore deleted afte analysis
#        pb <- txtProgressBar(1,length(mRNA))
#        for (locus in 1:length(mRNA)) {
#            current <- mRNA[locus,]
#            seqsmRNA[[locus]] <- genomemod[[as.character(seqnames(current))]][start(current):end(current)] # this might be rewritten more elegantly...
#            setTxtProgressBar(pb, locus)
#        }
#        close(pb)
                                        #        save(seqsmRNA,file="seqsmRNA.rdata")
        load("/home/sm934/siRNAold/siRNA/seqsmRNA.rdata")
        gccontmRNA <- rowSums(alphabetFrequency(seqsmRNA)[,c("C","G")]) / rowSums(alphabetFrequency(seqsmRNA)[,c("C","G","A","T")]) # mean 0.37
        save(genes,annogr,miRNAs,transposonsgff,irs, ncRNA,mRNA,gccontmRNA,file="annotation_tigr.Rdata")
    }

library(MKmisc)
assignCI <- function(cis, windows, comma) {
    intwin <- cbind(findInterval(cis[1,], windows), findInterval(cis[2,], windows))
    z <- intwin[,2] - intwin[,1]
    asswin <- rep(NA, nrow(intwin))
    if(any(z == 0, na.rm = TRUE))
       asswin[which(z == 0)] <- intwin[which(z == 0), 1]
    if(any(z == 1, na.rm = TRUE))
        asswin[which(z == 1)] <- intwin[cbind(which(z == 1), apply(abs(intwin[which(z == 1),] - (comma + 1)), 1, which.min))]
    classIDs <- names(windows)[asswin]

    classIDs <- ordered(as.factor(classIDs), levels = names(windows)[-length(windows)])
#    classIDs <- apply(cbind(c(-Inf, log2(divisions)), c(log2(divisions), Inf)), 1, paste, collapse = ",")[asswin]
                                        #    classIDs    
}


classCI <- function(x1, x2, probs, divisions, comma, plot = TRUE, main = "", col = c("red", "blue", "green")) {
    if(missing(probs))
        probs <- 1 / (1 + 1/ divisions)
    windows <- c(low = -Inf, probs, infinite = Inf)

    if(is.null(names(probs))) names(windows) <- apply(cbind(c(-Inf, round(probs,3)), round(c(probs, Inf), 3)), 1, paste, collapse = "-")
    
    ratioCI <- apply(cbind(x1, x2), 1, function(x) binomCI(x[1], sum(x), method = "modified jeffreys")$CI)    

                                        #    plot(density(log10(ratioCI[1,]), na.rm = TRUE)); abline(v = log10(probs))
#    plot(density(ratioCI[1,], na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))
#    plot(density(log10(ratioCI[2,]), na.rm = TRUE)); abline(v = log10(probs))
#    plot(density(ratioCI[2,], na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))

    classIDs <- assignCI(ratioCI, windows, comma)

    addSegments <- function(adjFunc = log2, nsegs = 100)
        {
            ymax = par("usr")[4]
            cols <- col[as.integer(classIDs[1:nsegs])]; cols[is.na(cols)] <- "black"
            nz <- (1:nsegs)[ratioCI[1,1:nsegs] != 0  | ratioCI[2,1:nsegs] != 1]
            plotCI <- ratioCI[,nz]
            cols <- cols[nz]
            segments(x0 = adjFunc(plotCI[1,]), y0 = 1:nsegs / nsegs * ymax, x1 = adjFunc(plotCI[2,]), col = cols)
        }
    
    plot(density(log2(colMeans(ratioCI[,!is.na(classIDs)])), na.rm = TRUE), main = main); abline(v = log2(probs), lty = 2, col = "red")
    addSegments()
    
    plot(density(colMeans(ratioCI[,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1), main = main); abline(v = (probs), lty = 2, col = "red")
    addSegments(adjFunc = I)

    p2r <- function(p) 1 / ((1 / p) - 1)
    plot(density(log2(p2r(colMeans(ratioCI[,!is.na(classIDs)])))), main = main); abline(v = log2(p2r(probs)), lty = 2, col = "red")
    addSegments(adjFunc = log2)
    
    plot(density(log10((ratioCI[1,!is.na(classIDs)])), na.rm = TRUE), main = main); abline(v = log10(probs), lty = 2, col = "red")
    addSegments(adjFunc = log10)
#    plot(density((ratioCI[1,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs), lty = 2, col = "red")
#    addSegments(adjFunc = I)
    
#    plot(density(log10((ratioCI[2,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs), lty = 2, col = "red")
#    addSegments(adjFunc = log10)
    
#    plot(density((ratioCI[2,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs), lty = 2, col = "red")
#    addSegments(adjFunc = I)
    
    classIDs
}



# counting biases - uses the alignment data object used to run the segmentation
countingBiases <- function(locAnn, cl) {
    load("/home/tjh48/Code/segmentMap_III/aD_noJunk.RData") #aD #alignmentData
    colnames(values(aD@alignments))[2] <- "multireads"

    #samples <- read.delim("Arabidopsis_samples.txt", as.is = TRUE)
    
    samples <- read.table("/projects/sebastian_mueller/locus_map/library_annotation98.csv",sep="\t",header=TRUE)
    
    wt <- samples$Plant.Expt.Type %in% c("WT","Col/Col","dpi10")
    
    aDwidths <- width(aD@alignments)
    aDnormal <- aD[aDwidths >= 21 & aDwidths <= 24,]
    aDnormalnomulti <- aDnormal[!duplicated(as.character(aDnormal@alignments$tag)),]
    
#    save(aDnormal,file="aDnormal.rdata")
#    load("aDnormal.rdata")
    aDwidths <- width(aDnormal@alignments)
    
    firstNuc <- substr(aDnormal@alignments$tag,1,1)
    firstNucnomulti <- substr(aDnormalnomulti@alignments$tag,1,1)
                                        #proportion of sRNAs with a given 5'nuc
    table(firstNucnomulti)/length(firstNucnomulti)
                                        #    A    C    G    T 
                                        # 0.38 0.14 0.19 0.28

    # find normal ratio of first base nucleotides
    expectedRatio <-  (tapply(rowSums(aDnormalnomulti@data),firstNucnomulti,sum)/sum(aDnormalnomulti@data))
                                        #    A    C    G    T 
                                        # 0.37 0.13 0.23 0.28

    # get counts in each locus
    aDA <- aDnormal[firstNuc=="A",]
    countsA <- getCounts(segments=locAnn,aD=aDA,cl=cl); rm(aDA)
    locAnn$countsA <- rowSums(countsA[,wt])
    aDT <- aDnormal[firstNuc=="T",]
    countsT <- getCounts(segments=locAnn,aD=aDT,cl=cl); rm(aDT)
    locAnn$countsT <- rowSums(countsT[,wt])
    aDC <- aDnormal[firstNuc=="C",]
    countsC <- getCounts(segments=locAnn,aD=aDC,cl=cl); rm(aDC)
    locAnn$countsC <- rowSums(countsC[,wt])
    aDG <- aDnormal[firstNuc=="G",]
    countsG <- getCounts(segments=locAnn,aD=aDG,cl=cl); rm(aDG)
    locAnn$countsG <- rowSums(countsG[,wt])
    
    
    firstBase <- cbind(locAnn$countsA, locAnn$countsC, locAnn$countsG, locAnn$countsT)

    # assumes binomial distribution and looks for significant variation for each locus
    z <- pbinom(firstBase,
                matrix(rowSums(firstBase), ncol = ncol(firstBase), nrow = nrow(firstBase)),
                prob = matrix(expectedRatio, ncol = ncol(firstBase), nrow = nrow(firstBase), byrow = TRUE), lower.tail = FALSE)
    z[rowSums(firstBase) == 0,] <- NA
    zq <- matrix(p.adjust(z, method = "BH"), ncol = ncol(firstBase))
    zq[rowSums(firstBase) == 0,] <- 1
    pred <- apply(do.call("cbind", lapply(1:4, function(ii) c("", c("A", "C", "G", "T")[ii])[as.integer(zq[,ii] < 0.01) + 1])), 1, paste, collapse = "")
    pred[pred == ""] <- NA
    
    locAnn$predominant_5prime_letter <- as.factor(pred)
    

                                        # counting numbers of 21s and 24s. For arabidopsis, 21s includes 21s and 22s, 24s includes 23s and 24s. Check with Bruno to see if this is sensible in Chlamy.
    aDs <- aD[width(aD@alignments)<21,]
    countss <- getCounts(segments=locAnn,aD=aDs,cl=cl); rm(aDs)
    countsswt <- countss[,wt]
    aD21 <- aD[width(aD@alignments)==21 | width(aD@alignments)==22,]
    counts21 <- getCounts(segments=locAnn,aD=aD21,cl=cl); rm(aD21)
    counts21wt <- counts21[,wt]
    aD24 <- aD[width(aD@alignments)==23 | width(aD@alignments)==24,]
    counts24 <- getCounts(segments=locAnn,aD=aD24,cl=cl); rm(aD24)
    counts24wt <- counts24[,wt]
    aDb <- aD[width(aD@alignments)>24,]
    countsb <- getCounts(segments=locAnn,aD=aDb,cl=cl); rm(aDb)
    countsbwt <- countsb[,wt]
        
                                        #sRNA size ratio
    locAnn$ratio21vs24 <- log2((rowSums(counts21wt))/(rowSums(counts24wt)))
    locAnn$counts21 <- rowSums(counts21wt)
    locAnn$counts24 <- rowSums(counts24wt)
    locAnn$ratio21vs24[(locAnn$counts21+locAnn$counts24)<=5] <- NaN


    # this function computes confidence intervals on the ratio of 21s and 24s for each locus, and then uses that to put each locus in a window. This function will produce a plot of the density of the mean (and the log of the mean) of the confidence intervals; choose the 'probs' values in such a way to split the modes of the density plots.
    locAnn$ratio21vs24Class <- classCI(x1 = locAnn$counts21, x2 = locAnn$counts24, probs = c(med = 10^-1.5, high = 10^-0.5), comma = 2, main = "21/24 ratio")

    # as above, but now looking at strand biases (i.e., whether most of the sRNAs are on the same strand, or evenly split across strands)
    countsminus <- getCounts(segments=locAnn,aD=aDnormal[strand(aDnormal@alignments)=="-",],cl=cl)
    countsplus <- getCounts(segments=locAnn,aD=aDnormal[strand(aDnormal@alignments)=="+",],cl=cl)
    countsall <- getCounts(segments=locAnn,aD=aDnormal, cl=cl)
    countsallwt <- countsall[,wt]
    countsminuswt <- countsminus[,wt]
    countspluswt <- countsplus[,wt]
    

                                        #strand bias
    locAnn$countsminus <- rowSums(countsminuswt)
    locAnn$countsplus <- rowSums(countspluswt)
    locAnn$ratio_strand <- log2((rowSums(countsminuswt))/(rowSums(countspluswt)))
    ##setting low expressed loci ratio to NaN due to low confidence
    locAnn$ratio_strand[(locAnn$countsminus+locAnn$countsplus)<=5] <- NaN
    locAnn$ratio_strand_abs <- abs(locAnn$ratio_strand)
    
                                        # confidence intervals on strands, as before
    
    locAnn$ratio_strand_class <- classCI(locAnn$countsplus, locAnn$countsminus,
                                         probs = c(0.2,0.35,0.65,0.8), comma = 2, main = "Strand bias", col = c("red", "blue", "green", "blue", "red"))
    levels(locAnn$ratio_strand_class) <- c("strong bias", "med bias", "no bias", "med bias", "strong bias")
    

        ##computing repetetivness for each loci (total reads div by multi match corrected)
    matches <- aDnormal@alignments$multireads
    countsnormalwt <- getCounts(segments=locAnn,aD=aDnormal[,wt],cl=cl)
    aDnormal@data <-  aDnormal@data/matches

    countsnormalwtnorm <- getCounts(segments=locAnn,aD=aDnormal[,wt],cl=cl)

    repCI <- apply(cbind(rowSums(countsnormalwtnorm), rowSums(countsnormalwt))[1:100,], 1, function(x) binomCI(x[1], x[2], method = "modified jeffreys")$CI)
    
    locAnn$repetitiveness <- 1-(rowSums(countsnormalwtnorm)/rowSums(countsnormalwt))
    locAnn$repetitivenessClass <- as.ordered(cut(locAnn$repetitiveness,quantile(locAnn$repetitiveness,probs=seq(0.25,1,0.25),na.rm=TRUE),include.lowest=TRUE, labels=rev(c("high","median","low"))))

    locAnn$repetitivenessClass <- classCI(rowSums(countsnormalwt) - rowSums(countsnormalwtnorm), rowSums(countsnormalwtnorm), probs = c(med = 0.3, high = 0.6), comma = 1, main = "Repetitiveness")

    locAnn
}


# this is an easy one. Might look at density map of (log) lengths of loci to check that the intervals (ignore intervals_fine) appear more or less sensible. 
sizeClasses <- function(locAnn, intervals = c(0, 50, 2000, Inf)) {
    intervals_fine <- c(14:30,32,35,40,50,60,70,80,90,100,120,140,160,180,200,224,250,275,300,350,400,450,500,600,700,800,1000,1200,1500,1700,2000,2500,3000,3500,4000,5000,Inf)

    #intervals <- c(0,30,100,200,400,2000,Inf)
    
    locAnn$size <- width(locAnn)
    locAnn$sizeclass <- as.ordered(cut(width(locAnn),intervals))
    locAnn$sizeclass_fine <- as.ordered(cut(width(locAnn),intervals_fine))
    locAnn
}

makeChipData <- function() {
                                        #chip_H3K9me2 <- readGappedAlignments("/home/sm934/seqdata/chipseq/H3K9me2_SRR037794_multi.bam") #GappedAlignments

#using readBAM for further processing 
    chrs = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5","mitochondria","chloroplast")
    chrlens <- c(30427671, 19698289, 23459830, 18585056, 26975502,154478,366924)
    
    chip_H3K9me2aD <- readBAM(c("H3K9me2_SRR037794_multi.bam","H3K27me1_SRR037788_multi.bam","H3_SRR037787_multi.bam","MRNA-seq_SRR066816_multi.bam","H3K27me3_SRR037881_multi.bam"),libnames=c("H3K9me2","H3K27me1","H3","MRNA-seq","H3K27me3"),dir="/home/sm934/seqdata/chipseq",repl=c(1:5),chrs=chrs,chrlens=chrlens, maxlen = 200)
#    chip_H3K9me2aD@alignments$matches <- chip_H3K9me2aD@alignments$multireads
    levels(chip_H3K9me2aD@replicates)<-c("H3K9me2","H3K27me1","H3","MRNA-seq","H3K27me3")
    save(chip_H3K9me2aD,file="chip_H3K9me2aD.rdata")
    
    chip_H3K4aD <- readBAM(c("H3K4me2_SRR037789_multi.bam","H3K4me3_SRR037791_multi.bam","H3K27me3_SRR037881_multi.bam","H3_SRR037787_multi.bam","MRNA-seq_SRR066816_multi.bam"),libnames=c("H3K4me2","H3K4me3","H3K27me3","H3unique","mRNAseq_unique"),dir="/home/sm934/seqdata/chipseq",repl=c(1:5),chrs=chrs,chrlens=chrlens, maxlen = 200)
#    chip_H3K4aD@alignments$matches <- chip_H3K4aD@alignments$multireads

    # ?some of these were supposed to be unique but these files are not available?
    
    levels(chip_H3K4aD@replicates)<-c("H3K4me2","H3K4me3","H3K27me3","H3unique","mRNAseq_unique")
    save(chip_H3K4aD,file="chip_H3K4aD.rdata")
}

chipAnn <- function(locAnn) {

    load("chip_H3K9me2aD.rdata")
    load("chip_H3K4aD.rdata")
    
                                        #bamFls <- list.files("/home/sm934/seqdata/chipseq/ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP002/SRP002100", "multi.bam$", full=TRUE)
    bamFls <- list.files("/home/sm934/seqdata/chipseq", "multi.bam$", full=TRUE)
    names(bamFls) <- sub("\\..*", "", basename(bamFls))
                                        #names(bamFlsunique) <- sub("\\..*", "", basename(bamFlsunique))
    counts_marks <- matrix(NA,nr=length(locAnn),ncol=length(bamFls))
    #counts_marks_unique <- matrix(NA,nr=length(locAnn),ncol=length(bamFlsunique))
    colnames(counts_marks) <- names(bamFls)
    #colnames(counts_marks_unique) <- names(bamFlsunique)
    rownames(counts_marks) <- names(locAnn)
                                        # rownames(counts_marks_unique) <- names(locAnn)
    
                                        #libsize_counts_marks <- rep(NA,length(bamFls))
                                        #names(libsize_counts_marks) <- names(bamFls)
    
    for (i in 1:length(bamFls)) {
	message(bamFls[i])
	aln <- readGAlignments(bamFls[i])
        seqlevels(aln) <- gsub("Chr", "", seqlevels(aln))
	hits <- countOverlaps(aln, locAnn)
	counts_marks[,names(bamFls[i])] <- countOverlaps(locAnn, aln[hits==1])
                                        #libsize_counts_marks[i] <- length(aln)
    }
                                        #values taken from bowtie status output ( reads with at least one reported alignment)
    libsizes_counts_marks_bowtie <- read.csv("/home/sm934/seqdata/chipseq/libsizes_2mm.csv",header=FALSE,row.names=1)
    
    
#    libsize_counts_marks_unique <- rep(NA,length(bamFlsunique))
#    names(libsize_counts_marks_unique) <- names(bamFlsunique)
    
#    for (i in 1:length(bamFlsunique)) {
#	message(bamFlsunique[i])
#	aln <- readGAlignments(bamFlsunique[i])
#	hits <- countOverlaps(aln, locAnn)
#	counts_marks_unique[,names(bamFlsunique[i])] <- countOverlaps(locAnn, aln[hits==1])
#	libsize_counts_marks_unique[i] <- length(aln)
#    }
#    rm(bamFls,bamFlsunique,hits,aln,i)
    
    counts_marks_rpm <- t(t(counts_marks)/libsizes_counts_marks_bowtie[,1])*1e6
    counts_marks_rpkm <- counts_marks_rpm/(width(locAnn@ranges)/1000)
    
    save(counts_marks,counts_marks_rpm,counts_marks_rpkm,file="histonmark_counts_loci4_2mm.rdata")
}


histoneAnnotate <- function(locAnn, file = "histonmark_counts_loci4_2mm.rdata")
    {
        libsizes_counts_marks_bowtie <- read.csv("/home/sm934/seqdata/chipseq/libsizes_2mm.csv",header=FALSE,row.names=1)
        load(file)
        
        locAnn$H3counts <- counts_marks[,"H3_SRR037787_multi"]
        ciH3 <- (apply(cbind(locAnn$H3counts, width(locAnn)), 1, function(z)
            poisson.test(z[1], z[2])$conf.int
                        ))
        locAnn$h3Class <- assignCI(ciH3, windows = c(none = -Inf, low = 0.01, med = 1, high = 10, infinite = Inf), comma = 0)
        plot(density(log10(apply(ciH3,2,mean))), main = "H3 accumulation rate per base")
        abline(v = log10(c(0.01, 1, 10)), lty = 2, col = "red")
        segments(x0 = log10(ciH3[1,1:100]), y0 = par("usr")[4] * 1:100 / 100, x1 = log10(ciH3[2,1:100]), col = c("red", "blue", "green", "purple")[as.integer(locAnn$h3Class[1:100])])
        savePlot("h3_CI.png")
        
        locAnn$H3K27me1_counts <- counts_marks[,"H3K27me1_SRR037788_multi"]
        locAnn$H3K27me3_counts <- counts_marks[,"H3K27me3_SRR037881_multi"]
        locAnn$H3K36me2_counts <- counts_marks[,"H3K36me2_SRR037785_multi"]
        locAnn$H3K36me3_counts <- counts_marks[,"H3K36me3_SRR037786_multi"]
        locAnn$H3K4me2_counts <- counts_marks[,"H3K4me2_SRR037790_multi"]
        locAnn$H3K4me3_counts <- counts_marks[,"H3K4me3_SRR037792_multi"]
        locAnn$H3K9Ac_counts <- counts_marks[,"H3K9Ac_SRR037879_multi"]
        locAnn$H3K9me2_counts <- counts_marks[,"H3K9me2_SRR037793_multi"]
        locAnn$H3K18Ac_counts <- counts_marks[,"SRR037710_H3K18Ac_multi"]               

        locAnn$hasExpression <- factor(counts_marks[,"MRNA-seq_SRR066817_multi"] > 0)

        locAnn$hasH3K9me2 <- classCI(x1 = locAnn$H3K9me2_counts, locAnn$H3counts, probs = c(med = 10^-1, high = 10^-0.3), comma = 0, main = "H3K9me2")
        savePlot("H3K9me2_CI.png")
        
        locAnn$hasH3K27me1 <- classCI(x1 = locAnn$H3K27me1_counts, locAnn$H3counts, probs = c(med = 10^-1, high = 0.6), comma = 0, main = "H3K27me1"); savePlot("H3K27me1_CI.png")
        locAnn$hasH3K27me3 <- classCI(x1 = locAnn$H3K27me3_counts, locAnn$H3counts, probs = c(med = 10^-2, high = 10^-1.3), comma = 0, main = "H3K27me3"); savePlot("H3K27me3_CI.png")
        locAnn$hasH3K36me2 <- classCI(x1 = locAnn$H3K36me2_counts, locAnn$H3counts, probs = c(med = 10^-1.7, high = 10^-0.6), comma = 0, main = "H3K36me2"); savePlot("H3K36me2_CI.png")
        locAnn$hasH3K36me3 <- classCI(x1 = locAnn$H3K36me3_counts, locAnn$H3counts, probs = c(med = 10^-2, high = 10^-0.8), comma = 0, main = "H3K36me3"); savePlot("H3K36me3_CI.png")
        locAnn$hasH3K4me2 <- classCI(x1 = locAnn$H3K4me2_counts, locAnn$H3counts, probs = c(med = 10^-2.8, high = 10^-1.2), comma = 0, main = "H3K4me2"); savePlot("H3K4me2_CI.png")
        locAnn$hasH3K4me3 <- classCI(x1 = locAnn$H3K4me3_counts, locAnn$H3counts, probs = c(med = 10^-2, high = 10^-1), comma = 0, main = "H3K4me3"); savePlot("H3K4me3_CI.png")
        locAnn$hasH3K9Ac <- classCI(x1 = locAnn$H3K9Ac_counts, locAnn$H3counts, probs = c(high = 10^-1.1), comma = 0, main = "H3K9Ac"); savePlot("H3K9Ac_CI.png")
        locAnn
    }


methAnnotate <- function(locAnn, cl)
    {
        library(stats4)

                
        ###########bringing in methylation marks (Stroud, Cell, '13)##############
        ##for each sequencing context, the proportions of methylated C is calculated ranging
        ##from 0 (no C in this context is methylated at all in this locus, to 1, all C are methylated all the time)
        
        meth_CGwt=import("/data/public_data/arabidopsis/Stroud/GSM980986_WT_rep2_CG.bwig", format = "BigWig")
        meth_CHHwt=import("/data/public_data/arabidopsis/Stroud/GSM980986_WT_rep2_CHH.bwig", format = "BigWig")
        meth_CHGwt=import("/data/public_data/arabidopsis/Stroud/GSM980986_WT_rep2_CHG.bwig", format = "BigWig")
        
        seqlevels(meth_CGwt) <- seqlevels(meth_CHGwt) <- seqlevels(meth_CHHwt) <- seqlevels(locAnn)
        
        tmp=findOverlaps(meth_CGwt,locAnn)
        meth_CGwtloci <- rep(NA,length(locAnn))
        methCGloc <- list()
        for (i in unique(subjectHits(tmp))) {
                                        #meth_CGwtloci[i] <- mean(abs(meth_CGwt[queryHits(tmp)[subjectHits(tmp)==i],]$score))
            methCGloc[[i]] <- abs(meth_CGwt[queryHits(tmp)[subjectHits(tmp)==i],]$score)
        }

        tmp=findOverlaps(meth_CHGwt,locAnn)
        meth_CHGwtloci <- rep(NA,length(locAnn))
        methCHGloc <- list()
        for (i in unique(subjectHits(tmp))) {
                                        #meth_CHGwtloci[i] <- mean(abs(meth_CHGwt[queryHits(tmp)[subjectHits(tmp)==i],]$score))
            methCHGloc[[i]] <- abs(meth_CHGwt[queryHits(tmp)[subjectHits(tmp)==i],]$score)
        }

        tmp=findOverlaps(meth_CHHwt,locAnn)
        meth_CHHwtloci <- rep(NA,length(locAnn))
        methCHHloc <- list()
        for (i in unique(subjectHits(tmp))) {
                                        #meth_CHHwtloci[i] <- mean(abs(meth_CHHwt[queryHits(tmp)[subjectHits(tmp)==i],]$score))
            methCHHloc[[i]] <- abs(meth_CHHwt[queryHits(tmp)[subjectHits(tmp)==i],]$score)
        }

        
#        methCG.ci <- sapply(methCGloc[1:450], function(x) {

#            nloglikbeta = function(mu, sig) {
#                alpha = mu^2*(1-mu)/sig^2-mu
#                beta = alpha*(1/mu-1)
#                -sum(dbeta(x, alpha, beta, log=TRUE))
#            }
#            
#            if(length(x) == 0) return(c(NA, NA))            
#            if(all(x == x[1])) return(c(x[1],x[1]))
#            x <- pmin(x, 0.99); x <- pmax(x, 0.01)
#            est = mle(nloglikbeta, start=list(mu=mean(x), sig=min(sd(x), (mean(x) - mean(x)^2) * 0.9)))
#            confint(est)[1,]
#        })

        bootMean <- function(x) {
            library(simpleboot)
            if(length(x) == 0) return(c(NA, NA))            
            if(all(x == x[1])) return(c(x[1],x[1]))            
            x.boot = one.boot(x, mean, R=10^4)
            boot.ci(x.boot, type="bca")$bca[,4:5]
        }
        
        methCG.bci <- parSapply(cl, methCGloc, bootMean)
        locAnn$CGclass <- assignCI(methCG.bci, windows = c(low = -Inf, med = 0.2, high = 0.6, inf = Inf), comma = 0)
        plot(density(apply(methCG.bci, 2, mean), na.rm = TRUE, from = 0, to = 1), main = "CpG context methylation"); abline(v = c(0.2, 0.6), lty = 2, col = "red")
        segments(x0 = methCG.bci[1,1:100], y0 = par("usr")[4] * 1:100 / 100, x1 = methCG.bci[2,1:100], col = c("red", "blue", "green")[as.integer(locAnn$CGclass[1:100])])
        savePlot("CGmet_CIs.png")
        
        methCHG.bci <- parSapply(cl, methCHGloc, bootMean)
        locAnn$CHGclass <- assignCI(methCHG.bci, windows = c(low = -Inf, high = 0.1, inf = Inf), comma = 0)
        plot(density(apply(methCHG.bci, 2, mean), na.rm = TRUE, from = 0, to = 1), main = "CHG context methylation", ); abline(v = c(0.1), lty = 2, col = "red")
        segments(x0 = methCHG.bci[1,1:100], y0 = par("usr")[4] * 1:100 / 100, x1 = methCHG.bci[2,1:100], col = c("red", "blue", "green")[as.integer(locAnn$CHGclass[1:100])])
        savePlot("CHGmet_CIs.png")

        methCHH.bci <- parSapply(cl, methCHHloc, bootMean)
        locAnn$CHHclass <- assignCI(methCHH.bci, windows = c(low = -Inf, high = 0.05, inf = Inf), comma = 0)
        plot(density(apply(methCHH.bci, 2, mean), na.rm = TRUE,from = 0, to = 1), main = "CHH context methylation", ); abline(v = 0.05, lty = 2, col = "red")
        segments(x0 = methCHH.bci[1,1:100], y0 = par("usr")[4] * 1:100 / 100, x1 = methCHH.bci[2,1:100], col = c("red", "blue", "green")[as.integer(locAnn$CHHclass[1:100])])
        savePlot("CHHmet_CIs.png")
        
        
        locAnn
}

                                        # This function needs an annotation map, here stored in 'annotation_tigr.Rdata'. This gets constructed in the TIGR function above

featureAnn <- function(locAnn, loci, file = "annotation_tigr.Rdata") {
    load(file)
##########genetic elements overlap##########
                                        #finding out about how often types overlap in general (if we get with ignoring amigously assigned loci to gene annotation) (seb_functions.R)
#    overlaptable <- cross_overlap_grange(genes)
                                        #most: Promotor / gene (8243), promotor / Helitron (791), Gene / Helitron (791), Gene / MudR (394)
    
                                        #prioritylist <- levels(genes$type)
                                        #reorder this variable manually using the following output as template
                                        #dump("prioritylist",file="") 
    prioritylist <-
        c("miRNA", "rRNA", "tRNA","NAT","ncRNA","lincRNA",  "gene",  "DNA/En-Spm", "DNA/Harbinger", "DNA/HAT", "DNA/Mariner", 
          "DNA/MuDR", "DNA/Pogo", "DNA/Tc1","DNA", "LINE?", "LINE/L1", "LTR/Copia", 
          "LTR/Gypsy", "RathE", "RC/Helitron", "SINE", "pseudogene","promotor","Unassigned")
    
    prioritylistte <-
        c("DNA/En-Spm", "DNA/Harbinger", "DNA/HAT", "DNA/Mariner", 
          "DNA/MuDR", "DNA/Pogo", "DNA/Tc1","DNA", "LINE?", "LINE/L1", "LTR/Copia", 
          "LTR/Gypsy", "RathE", "RC/Helitron", "SINE","Unassigned","miRNA", "rRNA", "tRNA","NAT","ncRNA","lincRNA",  "gene",   "pseudogene","promotor")
    
    try(rm(genesreordered))
    for (type in prioritylist) {if(exists("genesreordered")) {genesreordered <- c(genesreordered,genes[genes$type==type,]) } else {genesreordered <- genes[genes$type==type,]}}
    
    try(rm(genesreorderedte))
    for (type in prioritylistte) {if(exists("genesreorderedte")) {genesreorderedte <- c(genesreorderedte,genes[genes$type==type,]) } else {genesreorderedte <- genes[genes$type==type,]}}
    
    
#    genesreorderedcoarse <- genesreordered
#    levels(genesreorderedcoarse$type) <- c(levels(genesreorderedcoarse$type),"TE")
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="Unassigned"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/En-Spm"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/Harbinger"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/HAT"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/Mariner"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/MuDR"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/Pogo"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA/Tc1"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="DNA"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="LINE?"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="LINE/L1"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="LTR/Copia"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="LTR/Gypsy"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="RathE"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="RC/Helitron"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="SINE"] <- "TE"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="lincRNA"] <- "ncRNA"
#    genesreorderedcoarse$type[genesreorderedcoarse$type=="RC/Helitron"] <- "TE"
#    genesreorderedcoarse$type <- droplevels(genesreorderedcoarse$type)
    
#    load("/home/sm934/siRNAv4/genesTAIR10.rdata")
    
    table(genes$type)
                                        #overlap genes loci
    seqlevels(genesreordered)[1:5] <- seqlevels(locAnn)
    locivsgenes <- findOverlaps(locAnn,genesreordered)
    uniquegenes <- !duplicated(queryHits(locivsgenes))
    uniqueloci <- !duplicated(subjectHits(locivsgenes))


    # this bit finds the wild-type samples - insert your own
#    samples <- read.table("/projects/sebastian_mueller/locus_map/library_annotation98.csv",sep="\t",header=TRUE)
#    wt <- samples$Plant.Expt.Type %in% c("WT","Col/Col","dpi10")
                                        #tissue expression only for wt data
#    wtrepgroups <- as.integer(unique(loci@replicates[wt]))
    ##TODO: FDR?
#    locAnn$tissue_expression <- rowSums(loci@locLikelihoods[,wtrepgroups]>log(0.5))
#    locAnn$tissue_expressionClass <- as.ordered(cut(locAnn$tissue_expression,breaks=c(-Inf,1,10,Inf),include.lowest=TRUE,labels=c("specific","inbetween","common")))
    
    
    ##overlaptype vs loci

    ## this is the new bit - puts in individual columns for each annotation type.

    seqlevels(genes)[1:5] <- seqlevels(locAnn)
    for(lev in levels(values(genes)$type))
        values(locAnn)[,lev] <- getOverlaps(locAnn, genes[values(genes)$type == lev,], whichOverlaps = FALSE, cl = NULL)

                                        # Done. that was easy...
    
    locAnn$overlaptype <- rep("intergenic", length(locAnn))
    locAnn$overlaptype[queryHits(locivsgenes)[uniquegenes]] <- as.character(genesreordered[subjectHits(locivsgenes)[uniquegenes],]$type)
    locAnn$overlaptype[locAnn$overlaptype=="Unassigned"] <- "intergenic"
    
    seqlevels(transposonsgff) <- seqlevels(locAnn)
    teoverlap <- findOverlaps(locAnn,transposonsgff)
    teoverlapunique <- !rev(duplicated(rev(queryHits(teoverlap))))
    locAnn$isTE <- rep("none", length(locAnn))
    locAnn$isTE[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$type)[subjectHits(teoverlap)[teoverlapunique]]
    locAnn$overlaptype[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$type)[subjectHits(teoverlap)[teoverlapunique]]

                                        #inverted repeats overlap
    irs1000 <- irs[irs$score>1000,]
    seqlevels(irs1000) <- seqlevels(locAnn)
    iroverlap <- findOverlaps(locAnn,irs1000)
    iroverlapunique <- !rev(duplicated(rev(queryHits(iroverlap))))
    locAnn$isIR <- rep(FALSE, length(locAnn))
    locAnn$whichIR <- rep(FALSE, length(locAnn))
    locAnn$whichIR[queryHits(iroverlap)[iroverlapunique]] <- as.character(irs1000$ID)[subjectHits(iroverlap)[iroverlapunique]]
    locAnn$isIR[queryHits(iroverlap)[iroverlapunique]] <- TRUE
    locAnn$isIR <- ordered(locAnn$isIR)
    
                                        #reordering factor levels for overlaptype for better legends
    
    locAnn$overlaptype <- factor(locAnn$overlaptype, levels =
                                     c("intergenic","gene", "promotor", "ncRNA", "lincRNA","NAT", "miRNA", "rRNA", "tRNA", "pseudogene", "DNA", "DNA/En-Spm", "DNA/Harbinger", "DNA/HAT","DNA/Mariner", "DNA/MuDR", "DNA/Pogo", "DNA/Tc1", "RC/Helitron", "SINE", "LINE?", "LINE/L1", "LTR/Copia", "LTR/Gypsy", "RathE"))
    
    locAnn$overlaptypecoarse <- locAnn$overlaptype
    locAnn$overlaptypecoarse[locAnn$overlaptype=="lincRNA"] <- "ncRNA"
    locAnn$overlaptypecoarse[locAnn$overlaptype=="LINE?"] <- "LINE/L1"
    locAnn$overlaptypecoarse[locAnn$overlaptype=="DNA/Tc1"] <- "DNA"
    locAnn$overlaptypecoarse[locAnn$overlaptype=="tRNA"] <- "gene"
    locAnn$overlaptypecoarse[locAnn$overlaptype=="rRNA"] <- "gene"
    locAnn$overlaptypecoarse[locAnn$overlaptype=="pseudogene"] <- "ncRNA"
    locAnn$overlaptypecoarse <- droplevels(locAnn$overlaptypecoarse)

    locAnn
}

tissueSpec <- function(locAnn, loci) {
    selLoc <- selectLoci(loci, FDR = 0.1, perReplicate = TRUE, returnBool = TRUE)
                                        #AGOReps <- grep("AGO.*IP", samples[,2])
                                        #AGOReps <- unique(samples$Replicate.group[AGOReps])
                                        #AGOdat <- selLoc[,AGOReps]
                                        #
    samples <- read.table("/projects/sebastian_mueller/locus_map/library_annotation98.csv",sep="\t",header=TRUE)
    wt <- samples$Plant.Expt.Type %in% c("WT","Col/Col","dpi10")
                                        #tissue expression only for wt data
    wtrepgroups <- as.integer(unique(loci@replicates[wt]))
    
    ##tissue type expression
                                        #for (tissue in levels(classSegLike@annotation[wt,]$Tissue.type)[1:9]) {
    locAnn$"Arial" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Arial" & wt)]))]) > 0
    locAnn$"Floral" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Floral" & wt)]))]) > 0
    locAnn$"Floral/Silique" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Floral/Silique" & wt)])), drop = FALSE]) > 0
    locAnn$"Inflorescence" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Inflorescence" & wt)]))]) > 0
    locAnn$"Leaf" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Leaf" & wt)]))]) > 0
    locAnn$"Root" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Root" & wt)])),drop = FALSE]) > 0
    locAnn$"Seedling" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Seedling" & wt)]))]) > 0
    locAnn$"Shoot" <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Tissue.type == "Shoot" & wt)])),drop = FALSE]) > 0
    
                                        #floral specific
    locAnn$"Floralspecific" <- locAnn$"Floral" & !(locAnn$"Shoot" | locAnn$"Root" | locAnn$"Leaf" | locAnn$"Seedling"| locAnn$"Arial")
                                        #Gypsy/Copia/Tc1/En-Spm are enriched in Floral specific loci 
                                        #SINE,RathE,miRNA, rRNA is depleted
                                        #similar in Inflorescence!
                                        #Root
    locAnn$"Rootspecific" <- locAnn$"Root" & !(locAnn$"Shoot" | locAnn$"Floral"| locAnn$"Inflorescence"  | locAnn$"Leaf" | locAnn$"Seedling"| locAnn$"Arial")
    
                                        #hardly any there, few ncRNA/NAT bit specific
                                        #Arial
    
                                        #ncRNA/NAT/linc/DNA-Tc1/DNA-Pogo is specific
                                        #Leaf
    locAnn$"Leafspecific" <- locAnn$"Leaf" & !(locAnn$"Shoot" | locAnn$"Floral"| locAnn$"Inflorescence"  | locAnn$"Root" | locAnn$"Seedling"| locAnn$"Arial")
    locAnn
}


agoIP <- function(locAnn, loci) {
    ##assigning agoIPdata
    samples <- read.table("/projects/sebastian_mueller/locus_map/library_annotation98.csv",sep="\t",header=TRUE)
                                        #samples = classSegLike@annotation
    
                                        #hardly any, Gene sRNAs and NAT
    
                                        #tissue expression only for wt data    
    selLoc <- selectLoci(loci, FDR = 0.05, perReplicate = TRUE, returnBool = TRUE)

    
    ##incorperating Ago IP data:
    ##TODO, FDR cutoff
    #locAnn$Ago1IPclass <- as.factor(selLoc[,27])	# SL44/50, natasha, leaf
    #locAnn$Ago1IPbclass <- as.factor(selLoc[,25])	# SL35, susi, floral
    #locAnn$Ago1IPclass <- as.factor(selLoc[,43])	# SL206, susi, mock
    #locAnn$Ago1IPdclass <- as.factor(selLoc[,3])	# SL34, susi, florak

    ago1reps <- unique(samples[grep("AGO1 IP", samples[,2]), 5])
    ago1loc <- loci[,loci@replicates %in% ago1reps]; ago1loc@locLikelihoods <- ago1loc@locLikelihoods[,ago1reps]
    locAnn$Ago1IPMclass <-
        selectLoci(ago1loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)

    ago2reps <- unique(samples[grep("AGO2 IP", samples[,2]), 5])
    ago2loc <- loci[,loci@replicates %in% ago2reps]; ago2loc@locLikelihoods <- ago2loc@locLikelihoods[,ago2reps]
    locAnn$Ago2IPMclass <-
        selectLoci(ago2loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)

    ago4reps <- unique(samples[grep("AGO4 IP", samples[,2]), 5])
    ago4loc <- loci[,loci@replicates %in% ago4reps]; ago4loc@locLikelihoods <- ago4loc@locLikelihoods[,ago4reps]
    locAnn$Ago4IPMclass <-
        selectLoci(ago4loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)

    ago5reps <- unique(samples[grep("AGO5 IP", samples[,2]), 5])
    ago5loc <- loci[,loci@replicates %in% ago5reps]; ago5loc@locLikelihoods <- ago5loc@locLikelihoods[,ago5reps,drop = FALSE]
    locAnn$Ago5IPMclass <-
        selectLoci(ago5loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)

    ago6reps <- unique(samples[grep("AGO6 IP", samples[,2]), 5])
    ago6loc <- loci[,loci@replicates %in% ago6reps]; ago6loc@locLikelihoods <- ago6loc@locLikelihoods[,ago6reps]
    locAnn$Ago6IPMclass <-
        selectLoci(ago6loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)

    ago9reps <- unique(samples[grep("AGO9 IP", samples[,2]), 5])
    ago9loc <- loci[,loci@replicates %in% ago9reps]; ago9loc@locLikelihoods <- ago9loc@locLikelihoods[,ago9reps]
    locAnn$Ago9IPMclass <-
        selectLoci(ago9loc, FDR = 0.05, perReplicate = FALSE, returnBool = TRUE)
        
    locAnn$dclTripleClass <- as.factor(selLoc[,39])	# SL300, atilla, shoot
    
                                        #locAnn$Ago1IPclass <- (locAnn$Ago1IP + locAnn$Ago1IPb + locAnn$Ago1IPc + locAnn$Ago1IPd )/4 > log(0.7)
                                        #locAnn$Ago2IPclass <- ( locAnn$Ago2IP + locAnn$Ago2IPb )/2 > log(0.8)
                                        #locAnn$Ago4IPclass <- (locAnn$Ago4IP + locAnn$Ago4IPb + locAnn$Ago4IPc)/3 > log(0.7)
    locAnn
    
}

agoIP2 <- function(locAnn) {
    load("agoCD_pairwiseDE.RData")
    samples <- read.table("/projects/sebastian_mueller/locus_map/library_annotation98.csv",sep="\t",header=TRUE)
    agoID <- do.call("cbind.data.frame", lapply(agoCDLs, function(cdl) {
        agP <- rep(FALSE, length(locAnn))
        agP[match(selectTop(cdl, group = 2,
                  ordering = paste(unique(cdl@groups[[2]][which(samples[match(cdl@replicates, samples[,5]),2] != "WT")]), ">",
                      unique(cdl@groups[[2]][which(samples[match(cdl@replicates, samples[,5]),2] == "WT")]), sep = "")
                          , FDR = 0.05)@coordinates, locAnn)] <- TRUE
        as.factor(agP)
    })
                     )
    colnames(agoID) <- names(agoCDLs)
    values(locAnn) <- cbind(values(locAnn), agoID)
    locAnn
}
    
    
Pol45 <-  function(locAnn, cl) {
    ##defining RDR2 dependent loci: should be RdDM loci!
                                        #importing mutants from Lee12
    libfilesLee <- c(
        "Lee12_GSM893112_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893113_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893114_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893115_nrpe_1_adapters_trimmed_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893116_nrpe_2_adapters_trimmed_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893117_nrpe_3_adapters_trimmed_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893124_rdr2_adapters_trimmed_Inflorescence_redundant_multi.bam",
        "Lee12_GSM893123_nrpd1_adapters_trimmed_Inflorescence_redundant_multi.bam")

    chrs = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5","mitochondria","chloroplast")
    chrlens <- c(30427671, 19698289, 23459830, 18585056, 26975502,154478,366924)
    
    libnames <- paste(c(rep("WT",3),rep("nrpe",3),"rdr2","nrpd1"),c(1:3,1:3,1,1),sep="")
                                        #important: XS flag for redundant reads
    readsLee <- readBAM(libfilesLee,libnames=libnames,dir="/home/sm934/seqdata/sRNAexternallibs",repl=c(rep(c("WT"),3),rep(c("nrpe"),3),c("rdr2","nrpd1")),chrs=chrs,chrlens=chrlens)
                                        #load("/home/sm934/siRNAv2/readsLee12.rdata")
    seqlevels(readsLee@alignments)[1:5] <- seqlevels(locAnn)
    countsLee <- getCounts(segments=locAnn,aD=readsLee,cl=cl)
    colnames(countsLee) <- colnames(readsLee@data)
    
    ##testing PolV
    CDloci <- new("countData", data = countsLee[,1:6], replicates = c(rep(1,3), rep(2,3)), groups =  list(NDE = c(1,1,1,1,1,1), DE = c(rep(1,3), rep(2,3))),libsizes=readsLee@libsizes[1:6],annotation=as.data.frame(locAnn))
    
    CDlociP <- getPriors.NB(CDloci, samplesize = 1e5, estimation = "QL", cl = cl)
    CDlociP <- makeOrderings(CDlociP)
    CDlociL <- getLikelihoods(CDlociP, pET = 'BIC', cl = cl, nullData = TRUE, modelPriorSets = split(1:nrow(CDlociP), CDlociP@orderings[,2]))
    
                                        #tcloci <- topCounts(CDlociL, group = "DE",number=length(locAnn),normaliseData = TRUE)    
    
    ##extracting RdDM loci altogether (Pol4 and RDR2 mutants are pooled here to check loci for absence of smallRNAs compared to WT) 

    CDlocipol4 <- new("countData", data = countsLee[,c(1:3,8)], replicates = c(rep(1,3), 2), groups =  list(NDE = c(1,1,1,1), DE = c(rep(1,3), rep(2,1))),libsizes=readsLee@libsizes[c(1:3,8)],annotation=as.data.frame(locAnn))
        
    CDlociPpol4 <- getPriors.NB(CDlocipol4, samplesize = 1e5, estimation = "QL", cl = cl)
    CDlociPpol4 <- makeOrderings(CDlociPpol4)
    CDlociLpol4 <- getLikelihoods(CDlociPpol4, pET = 'BIC', cl = cl, nullData = TRUE, modelPriorSets = split(1:nrow(CDlociPpol4), CDlociPpol4@orderings[,2]))

    CDlocirdr2 <- new("countData", data = countsLee[,c(1:3,7)], replicates = c(rep(1,3), 2), groups =  list(NDE = c(1,1,1,1), DE = c(rep(1,3), rep(2,1))),libsizes=readsLee@libsizes[c(1:3,7)],annotation=as.data.frame(locAnn))
        
    CDlociPrdr2 <- getPriors.NB(CDlocirdr2, samplesize = 1e5, estimation = "QL", cl = cl)
    CDlociPrdr2 <- makeOrderings(CDlociPrdr2)
    CDlociLrdr2 <- getLikelihoods(CDlociPrdr2, pET = 'BIC', cl = cl, nullData = TRUE, modelPriorSets = split(1:nrow(CDlociPrdr2), CDlociPrdr2@orderings[,2]))

    
    save(countsLee,readsLee,CDlociL,CDlociLpol4,CDlociLrdr2, file="readsLee12.rdata")
}


annPol <- function(locAnn)
    {
        load("readsLee12.rdata")
        locAnn$polV_new_DE_posterior <- CDlociL@posteriors[,"DE"]
        locAnn$polV_new <- rep("not_known",length(locAnn)) 
        locAnn$polV_new[rownames(CDlociL@annotation) %in% rownames(selectTop(CDlociL, group = "DE", ordering = "1>2", FDR = 0.05)@annotation)] <- "PolV_dependent"
        pg <- exp(CDlociL@posteriors[,2]) * as.integer(CDlociL@orderings[,2] != "1>2")
        CDlociL@posteriors[,1] <- log(exp(CDlociL@posteriors[,1]) + pg)
        CDlociL@posteriors[,2] <- log(exp(CDlociL@posteriors[,2]) - pg)
        locAnn$polV_new[rownames(CDlociL@annotation) %in% rownames(selectTop(CDlociL, group = "NDE", FDR = 0.05)@annotation)] <- "PolV_independent"
        locAnn$polV_new <- ordered(locAnn$polV_new,levels = c("PolV_independent","not_known","PolV_dependent"))
#        tapply(locAnn$overlaptype,locAnn$polV,table)
                                        #RC/Helitron,RathE*,  LINE/L1,DNA/HAT,DNA/Harbinger  mostly polV dependent
                                        # a bit promoters, (like in zhong 12)
                                        #Gypsy mostly polV independent (see also zemach 13) also DNA/En-Spm,DNA,tRNA
        
                                        #68% estimated to be DE, so p-values shoud be interpreted carefully
                                        #upregulated loci in RdDM mutants are strand specific, have a diverse range of smallRNA width, and tend to associate with rRNA        

        locAnn$polIV_new <- rep("not_known",length(locAnn)) 
        locAnn$polIV_new[rownames(CDlociLpol4@annotation) %in% rownames(selectTop(CDlociLpol4, group = "DE", ordering = "1>2", FDR = 0.05)@annotation)] <- "PolIV_dependent"
        pg <- exp(CDlociLpol4@posteriors[,2]) * as.integer(CDlociLpol4@orderings[,2] != "1>2")
        CDlociLpol4@posteriors[,1] <- log(exp(CDlociLpol4@posteriors[,1]) + pg)
        CDlociLpol4@posteriors[,2] <- log(exp(CDlociLpol4@posteriors[,2]) - pg)
        locAnn$polIV_new[rownames(CDlociLpol4@annotation) %in% rownames(selectTop(CDlociLpol4, group = "NDE", FDR = 0.05)@annotation)] <- "PolIV_independent"
        locAnn$polIV_new <- ordered(locAnn$polIV_new,levels = c("PolIV_independent","not_known","PolIV_dependent"))

        
        locAnn$rdr2_new <- rep("not_known",length(locAnn)) 
        locAnn$rdr2_new[rownames(CDlociLrdr2@annotation) %in% rownames(selectTop(CDlociLrdr2, group = "DE", ordering = "1>2", FDR = 0.05)@annotation)] <- "Rdr2_dependent"
        pg <- exp(CDlociLrdr2@posteriors[,2]) * as.integer(CDlociLrdr2@orderings[,2] != "1>2")
        CDlociLrdr2@posteriors[,1] <- log(exp(CDlociLrdr2@posteriors[,1]) + pg)
        CDlociLrdr2@posteriors[,2] <- log(exp(CDlociLrdr2@posteriors[,2]) - pg)
        locAnn$rdr2_new[rownames(CDlociLrdr2@annotation) %in% rownames(selectTop(CDlociLrdr2, group = "NDE", FDR = 0.05)@annotation)] <- "Rdr2_independent"
        locAnn$rdr2_new <- ordered(locAnn$rdr2_new,levels = c("Rdr2_independent","not_known","Rdr2_dependent"))
        
                                        #as expect, almost no PolV dependent but RdDM indep loci!
        locAnn
    }
