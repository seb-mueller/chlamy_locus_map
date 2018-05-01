#Functions for use in annotation
#Adapted by Nick Matthews from code written by Tom Hardcastle and Sebastian Muller
#Date: 09/02/16

#Divides up loci into seperate size classes - use loci size plot to determain appropriate divisions
sizeClasses <- function(locAnn) {
      
    intervals <- c(0,30,75,150,1000,Inf) #different to 
    
    locAnn$size <- width(locAnn)
    locAnn$sizeclass <- as.ordered(cut(width(locAnn),intervals))
    locAnn
}

#Processes and generate annotations into single file
TIGR <- function()
    {
        transposonsgff <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/transposons_newcut.gff3") #GRanges
        SoloLTR <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/SoloLTR.gff3")
                                        #lincRNAs
        #lincRNAs <- import.gff3("/home/sm934/seqdata/lincRNAs_Liu/LincRNAs_identified_by_both_RepTAS_and_RNA-seq.gff")
        #lincRNAs$Name <- lincRNAs$ID
        #levels(lincRNAs$type) <- "lincRNA"
        
                                        #natRNAs
        #natRNAs <- import.gff3("/home/sm934/seqdata/Luo2013_NATs_stringentb.gff3")[,1:6]
        #natRNAs$type <- "NAT"
        

        irs <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/irf_output_Creinhardtii_236.gff")
	irs$type<-irs$source
	trs <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/trf_output_Creinhardtii_236.gff")
	trs$type<-trs$source       
                                        
        anno <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/Creinhardtii_281_v5.5.gene_exons.gff3")
                                        #annogr <- GRanges(seqnames=anno$space,ranges = IRanges(start=start(anno),end=end(anno)),strand=strand(anno),type=anno$type,group=anno$group)
        #annogr <- anno[!anno@elementMetadata$type=="chromosome",]
	#annogr<-anno
        
        miRNA <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/mirna_final_coordinates.gff")
	miRNA$type <- "miRNA"
        rRNA <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/Creinhardtii_rRNA.gff3")
	MSAT <- import.gff3("/home/bioinf/nem34/segmentMap_II/Old_Annotation_Files/Creinhardtii_MSAT.gff3")   
        #ncRNA <- annogr[annogr@elementMetadata$type=="ncRNA",] #align by blast - just same as transposons
        mRNA <- anno[anno$type=="mRNA",]
        genes<-anno[anno$type=="gene",]
	exons<-anno[anno$type=="exon",]
	CDS<-anno[anno$type=="CDS",]
	threeprimeUTR<-anno[anno$type=="three_prime_UTR",]
	fiveprimeUTR<-anno[anno$type=="five_prime_UTR",]
	
	
#I've made genes just the genes	- previous slightly unusual gene definition
	#elements <- annogr[annogr@elementMetadata$type=="miRNA" |annogr@elementMetadata$type=="ncRNA" | annogr@elementMetadata$type=="gene" | annogr@elementMetadata$type=="pseudogene" | annogr@elementMetadata$type=="rRNA" |  annogr@elementMetadata$type=="tRNA",] #only 11 sno,sn RNA, therefore discarded
        #elements$Index <- NULL;elements$Note <- NULL;elements$Parent <- NULL;elements$Index <- NULL;elements$Derives_from <- NULL;elements$Alias <- NULL
        #tRNA <- annogr[annogr@elementMetadata$type=="tRNA",] #don't have
        
     #Assume promoter region the 500bp flanking region around gene   
        promoter <- flank(anno[anno$type=="gene"],500)
        promoter$type <- droplevels(promoter$type)
        levels(promoter$type) <- "promoter"
#I've made genes just the genes - previous slightly unusual gene definition
	#genes <- c(elements,transposonsgff,lincRNAs,promoter,natRNAs)
        #genes$type <- droplevels(genes$type)
        #genes$sizeclass <- cut(width(genes),breaks=c(0,500,2000,max(width(genes))))
        #levels(genes$sizeclass) <- c("Short","Medium","Long")
        #genes$isTE <- genes$type %in% levels(transposonsgff$type) #not sure if will work or relevent - none fit at the moment so not using
	everything<-c(irs[,1:2],trs[,1:2],miRNA[,1:2],rRNA[,1:2],MSAT[,1:2],mRNA[,1:2],genes[,1:2],exons[,1:2],CDS[,1:2],threeprimeUTR[,1:2],fiveprimeUTR[,1:2],promoter[,1:2],transposonsgff[,1:2],SoloLTR[,1:2])
        
	save(genes,anno,transposonsgff,SoloLTR,rRNA,miRNA,MSAT,mRNA,irs,trs,promoter,threeprimeUTR,fiveprimeUTR,exons,CDS,everything,file="annotation_tigr.Rdata")
    }


#Function to comute overlaps with annotations using annotation file created by TIGR function

featureAnn <- function(locAnn, loci) {
    load("annotation_tigr.Rdata")
##########genetic elements overlap##########
                                        #finding out about how often types overlap in general (if we get with ignoring amigously assigned loci to gene annotation) (seb_functions.R)
    
        
    table(genes$type)
                                        #overlap genes loci
    geneoverlap <- findOverlaps(locAnn,genes) #got rid of genes reordered
    geneoverlapunique <- !duplicated(queryHits(geneoverlap))
    #uniqueloci <- !duplicated(subjectHits(locivsgenes))
##overlaptype vs loci
    #locAnn$overlaptype <- rep("intergenic", length(locAnn))
    #locAnn$overlaptype[queryHits(geneoverlap)[geneoverlapunique]] <- as.character(genes[subjectHits(geneoverlap)[geneoverlapunique],]$type)
    #locAnn$overlaptype[locAnn$overlaptype=="Unassigned"] <- "intergenic"
    locAnn$gene <- rep(FALSE, length(locAnn))
    locAnn$gene[queryHits(geneoverlap)[geneoverlapunique]] <- TRUE

  #seqlevels(genes)[1:5] <- seqlevels(locAnn) ##force same chr name - not 1:5 as that's for arabidopsis
  #  for(lev in levels(values(everything)$type)) 
  #      values(locAnn)[,lev] <- getOverlaps(locAnn, genes[values(genes)$type == lev,], whichOverlaps = FALSE, cl = NULL)
	

    # this bit finds the wild-type samples
    samples <- read.csv("/home/bioinf/nem34/segmentMap_II/Summary_of_Data_5.csv",header=TRUE)
    wt <- samples$GenuineControls %in% c("wt")
                                        #expression only for wt data
    wtrepgroups <- as.integer(unique(loci@replicates[wt]))
    locAnn$expression <- rowSums(loci@locLikelihoods[,wtrepgroups]>log(0.5))
    locAnn$expressionClass <- as.ordered(cut(locAnn$expression,breaks=c(-Inf,1,5,Inf),include.lowest=TRUE,labels=c("specific","inbetween","common"))) #changed common threshold to more than 5
    
        
    seqlevels(transposonsgff) <- seqlevels(locAnn)
    teoverlap <- findOverlaps(locAnn,transposonsgff)
    teoverlapunique <- !rev(duplicated(rev(queryHits(teoverlap))))
    locAnn$TE <- rep("none", length(locAnn))
    locAnn$TE[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$type)[subjectHits(teoverlap)[teoverlapunique]]
    #locAnn$overlaptype[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$type)[subjectHits(teoverlap)[teoverlapunique]]
    
    teoverlapunique <- !rev(duplicated(rev(queryHits(teoverlap))))
    locAnn$TEclass <- rep("none", length(locAnn))
    locAnn$TEclass[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$class)[subjectHits(teoverlap)[teoverlapunique]]
    #locAnn$overlaptype[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$class)[subjectHits(teoverlap)[teoverlapunique]]

	    teoverlapunique <- !rev(duplicated(rev(queryHits(teoverlap))))
    locAnn$TEorder <- rep("none", length(locAnn))
    locAnn$TEorder[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$order)[subjectHits(teoverlap)[teoverlapunique]]
    #locAnn$overlaptype[queryHits(teoverlap)[teoverlapunique]] <- as.character(transposonsgff$order)[subjectHits(teoverlap)[teoverlapunique]]
	
                                        #inverted repeats overlap - need to make sure I have
    seqlevels(irs) <- seqlevels(locAnn)
    iroverlap <- findOverlaps(locAnn,irs)
    iroverlapunique <- !rev(duplicated(rev(queryHits(iroverlap))))
    locAnn$IR <- rep(FALSE, length(locAnn))
    locAnn$IR[queryHits(iroverlap)[iroverlapunique]] <- TRUE


    seqlevels(trs) <- seqlevels(locAnn)
    troverlap <- findOverlaps(locAnn,trs)
    troverlapunique <- !rev(duplicated(rev(queryHits(troverlap))))
    locAnn$TR <- rep(FALSE, length(locAnn))
    locAnn$TR[queryHits(troverlap)[troverlapunique]] <- TRUE

    
#need to include miRNA,rRNAs, MSAT 
    seqlevels(MSAT) <- seqlevels(locAnn)
    MSAToverlap <- findOverlaps(locAnn,MSAT)
    MSAToverlapunique <- !rev(duplicated(rev(queryHits(MSAToverlap))))
    locAnn$MSAT <- rep(FALSE, length(locAnn))
    locAnn$MSAT[queryHits(MSAToverlap)[MSAToverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(MSAToverlap)[MSAToverlapunique]] <- as.character(MSAT$type)[subjectHits(MSAToverlap)[MSAToverlapunique]]	
		#MSATs only 2 hits...

    seqlevels(rRNA) <- seqlevels(locAnn)
    rRNAoverlap <- findOverlaps(locAnn,rRNA)
    rRNAoverlapunique <- !rev(duplicated(rev(queryHits(rRNAoverlap))))
    locAnn$rRNA <- rep(FALSE, length(locAnn))
    locAnn$rRNA[queryHits(rRNAoverlap)[rRNAoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(rRNAoverlap)[rRNAoverlapunique]] <- as.character(rRNA$type)[subjectHits(rRNAoverlap)[rRNAoverlapunique]]

    seqlevels(miRNA) <- seqlevels(locAnn)
    miRNAoverlap <- findOverlaps(locAnn,miRNA)
    miRNAoverlapunique <- !rev(duplicated(rev(queryHits(miRNAoverlap))))
    locAnn$miRNA <- rep(FALSE, length(locAnn))
    locAnn$miRNA[queryHits(miRNAoverlap)[miRNAoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(miRNAoverlap)[miRNAoverlapunique]] <- as.character(miRNA$type)[subjectHits(miRNAoverlap)[miRNAoverlapunique]]

#including promoters, UTRs, CDS, exons, introns...
    seqlevels(promoter) <- seqlevels(locAnn)
    promoteroverlap <- findOverlaps(locAnn,promoter)
    promoteroverlapunique <- !rev(duplicated(rev(queryHits(promoteroverlap))))
    locAnn$promoter <- rep(FALSE, length(locAnn))
    locAnn$promoter[queryHits(promoteroverlap)[promoteroverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(promoteroverlap)[promoteroverlapunique]] <- as.character(promoter$type)[subjectHits(promoteroverlap)[promoteroverlapunique]]

    seqlevels(mRNA) <- seqlevels(locAnn)
    mRNAoverlap <- findOverlaps(locAnn,mRNA)
    mRNAoverlapunique <- !rev(duplicated(rev(queryHits(mRNAoverlap))))
    locAnn$mRNA <- rep(FALSE, length(locAnn))
    locAnn$mRNA[queryHits(mRNAoverlap)[mRNAoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(mRNAoverlap)[mRNAoverlapunique]] <- as.character(mRNA$type)[subjectHits(mRNAoverlap)[mRNAoverlapunique]]

   
    seqlevels(CDS) <- seqlevels(locAnn)
    CDSoverlap <- findOverlaps(locAnn,CDS)
    CDSoverlapunique <- !rev(duplicated(rev(queryHits(CDSoverlap))))
    locAnn$CDS <- rep(FALSE, length(locAnn))
    locAnn$CDS[queryHits(CDSoverlap)[CDSoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(CDSoverlap)[CDSoverlapunique]] <- as.character(CDS$type)[subjectHits(CDSoverlap)[CDSoverlapunique]

    seqlevels(exons) <- seqlevels(locAnn)
    exonoverlap <- findOverlaps(locAnn,exons)
    exonoverlapunique <- !rev(duplicated(rev(queryHits(exonoverlap))))
    locAnn$exon <- rep(FALSE, length(locAnn))
    locAnn$exon[queryHits(exonoverlap)[exonoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(exonoverlap)[exonoverlapunique]] <- as.character(exon$type)[subjectHits(exonoverlap)[exonoverlapunique]

    seqlevels(fiveprimeUTR) <- seqlevels(locAnn)
    fiveprimeUTRoverlap <- findOverlaps(locAnn,fiveprimeUTR)
    fiveprimeUTRoverlapunique <- !rev(duplicated(rev(queryHits(fiveprimeUTRoverlap))))
    locAnn$fiveprimeUTR <- rep(FALSE, length(locAnn))
    locAnn$fiveprimeUTR[queryHits(fiveprimeUTRoverlap)[fiveprimeUTRoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(fiveprimeUTRoverlap)[fiveprimeUTRoverlapunique]] <- as.character(fiveprimeUTR$type)[subjectHits(fiveprimeUTRoverlap)[fiveprimeUTRoverlapunique]

    seqlevels(threeprimeUTR) <- seqlevels(locAnn)
    threeprimeUTRoverlap <- findOverlaps(locAnn,threeprimeUTR)
    threeprimeUTRoverlapunique <- !rev(duplicated(rev(queryHits(threeprimeUTRoverlap))))
    locAnn$threeprimeUTR <- rep(FALSE, length(locAnn))
    locAnn$threeprimeUTR[queryHits(threeprimeUTRoverlap)[threeprimeUTRoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(threeprimeUTRoverlap)[threeprimeUTRoverlapunique]] <- as.character(threeprimeUTR$type)[subjectHits(threeprimeUTRoverlap)[threeprimeUTRoverlapunique]

	#If want to include intron data:
    #seqlevels(introns) <- seqlevels(locAnn)
    #intronoverlap <- findOverlaps(locAnn,introns)
    #intronoverlapunique <- !rev(duplicated(rev(queryHits(intronoverlap))))
    #locAnn$isintron <- rep(FALSE, length(locAnn))
    #locAnn$isintron[queryHits(intronoverlap)[intronoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(intronoverlap)[intronoverlapunique]] <- as.character(intron$type)[subjectHits(intronoverlap)[intronoverlapunique]

#SoloLTRs	
	seqlevels(SoloLTR) <- seqlevels(locAnn)
    SoloLTRoverlap <- findOverlaps(locAnn,SoloLTR)
    SoloLTRoverlapunique <- !rev(duplicated(rev(queryHits(SoloLTRoverlap))))
    locAnn$SoloLTR <- rep(FALSE, length(locAnn))
    locAnn$SoloLTR[queryHits(SoloLTRoverlap)[SoloLTRoverlapunique]] <- TRUE
    #locAnn$overlaptype[queryHits(SoloLTRoverlap)[SoloLTRoverlapunique]] <- as.character(SoloLTR$type)[subjectHits(SoloLTRoverlap)[SoloLTRoverlapunique]

                                        #reordering factor levels for overlaptype for better legends
    #locAnn$overlaptype <- factor(locAnn$overlaptype, levels = c("intergenic","gene", "promoter", "Copia", "TOC1", "DIRS", "DNA", Gulliver", "Gypsy", "L1", "LTR", "Novosib", "P", "REM1", "REP", "RTE", "SINE", "TCR1", "TE", "TOC2"))
    
    locAnn
}

#Methylation Ovelaps - Meth IP data
methylation <- function(locAnn) {
    meth<-import.gff3("/data/pipeline/prod/SL55/SL55.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff3")
    seqlevels(meth) <- seqlevels(locAnn)
    methoverlap <- findOverlaps(locAnn,meth)
    methoverlapunique <- !rev(duplicated(rev(queryHits(methoverlap))))
    locAnn$meth <- rep(FALSE, length(locAnn))
    locAnn$meth[queryHits(methoverlap)[methoverlapunique]] <- TRUE
    locAnn$meth <- ordered(locAnn$meth)
    #locAnn$overlaptype[queryHits(methoverlap)[methoverlapunique]] <- as.character(meth$type)[subjectHits(methoverlap)[methoverlapunique]]

locAnn
}


#Methylation overlaps using new methylation loci from Tom
methAnnotate <- function(locAnn, cl)
    {
        library(stats4)
        
        methCG=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CGmeth.gff3")
        methCHH=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHHmeth.gff3")
        methCHG=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/chlamy_CHGmeth.gff3")
        
		seqlevels(methCG) <- seqlevels(methCHG) <- seqlevels(methCHH) <- seqlevels(locAnn)
        
		methCGoverlap <- findOverlaps(locAnn,methCG)
		methCGoverlapunique <- !rev(duplicated(rev(queryHits(methCGoverlap))))
		locAnn$methCG <- rep(FALSE, length(locAnn))
		locAnn$methCG[queryHits(methCGoverlap)[methCGoverlapunique]] <- TRUE
		locAnn$methCG <- ordered(locAnn$methCG)
		
	    methCHHoverlap <- findOverlaps(locAnn,methCHH)
		methCHHoverlapunique <- !rev(duplicated(rev(queryHits(methCHHoverlap))))
		locAnn$methCHH <- rep(FALSE, length(locAnn))
		locAnn$methCHH[queryHits(methCHHoverlap)[methCHHoverlapunique]] <- TRUE
		locAnn$methCHH <- ordered(locAnn$methCHH)
		
		methCHGoverlap <- findOverlaps(locAnn,methCHG)
		methCHGoverlapunique <- !rev(duplicated(rev(queryHits(methCHGoverlap))))
		locAnn$methCHG <- rep(FALSE, length(locAnn))
		locAnn$methCHG[queryHits(methCHGoverlap)[methCHGoverlapunique]] <- TRUE
		locAnn$methCHG <- ordered(locAnn$methCHG)
 
        locAnn
}


#Meth overlaps using differential loci from Tom - too small datasets
methDiff <- function(locAnn, cl)
    {
        library(stats4)
        
        methCG=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/Andy_differentialMeth_CG.gff3")
        methCHH=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/Andy_differentialMeth_CHH.gff3")
        methCHG=import.gff3("/home/bioinf/nem34/segmentMap_II/meth_data/Andy_differentialMeth_CHG.gff3")
        
		seqlevels(methCG) <- seqlevels(methCHG) <- seqlevels(methCHH) <- seqlevels(locAnn)
        
		#for CG meth
		methCG_wtonly <- methCG[methCG$ordering == "2>1"]
		methCG_dcl3only <- methCG[methCG$ordering == "1>2"]
		
		methCG_wtonlyoverlap <- findOverlaps(locAnn,methCG_wtonly)
		methCG_wtonlyoverlapunique <- !rev(duplicated(rev(queryHits(methCG_wtonlyoverlap))))
		locAnn$methCG_wtonly <- rep(FALSE, length(locAnn))
		locAnn$methCG_wtonly[queryHits(methCG_wtonlyoverlap)[methCG_wtonlyoverlapunique]] <- TRUE
		locAnn$methCG_wtonly <- ordered(locAnn$methCG_wtonly)
		
		methCG_dcl3onlyoverlap <- findOverlaps(locAnn,methCG_dcl3only)
		methCG_dcl3onlyoverlapunique <- !rev(duplicated(rev(queryHits(methCG_dcl3onlyoverlap))))
		locAnn$methCG_dcl3only <- rep(FALSE, length(locAnn))
		locAnn$methCG_dcl3only[queryHits(methCG_dcl3onlyoverlap)[methCG_dcl3onlyoverlapunique]] <- TRUE
		locAnn$methCG_dcl3only <- ordered(locAnn$methCG_dcl3only)
		
		#for CHH meth
		methCHH_wtonly <- methCHH[methCHH$ordering == "2>1"]
		methCHH_dcl3only <- methCHH[methCHH$ordering == "1>2"]
		
	    methCHH_wtonlyoverlap <- findOverlaps(locAnn,methCHH_wtonly)
		methCHH_wtonlyoverlapunique <- !rev(duplicated(rev(queryHits(methCHH_wtonlyoverlap))))
		locAnn$methCHH_wtonly <- rep(FALSE, length(locAnn))
		locAnn$methCHH_wtonly[queryHits(methCHH_wtonlyoverlap)[methCHH_wtonlyoverlapunique]] <- TRUE
		locAnn$methCHH_wtonly <- ordered(locAnn$methCHH_wtonly)
		
		methCHH_dcl3onlyoverlap <- findOverlaps(locAnn,methCHH_dcl3only)
		methCHH_dcl3onlyoverlapunique <- !rev(duplicated(rev(queryHits(methCHH_dcl3onlyoverlap))))
		locAnn$methCHH_dcl3only <- rep(FALSE, length(locAnn))
		locAnn$methCHH_dcl3only[queryHits(methCHH_dcl3onlyoverlap)[methCHH_dcl3onlyoverlapunique]] <- TRUE
		locAnn$methCHH_dcl3only <- ordered(locAnn$methCHH_dcl3only)
		
		#for CHG meth
		methCHG_wtonly <- methCHG[methCHG$ordering == "2>1"]
		methCHG_dcl3only <- methCHG[methCHG$ordering == "1>2"]
		
		methCHG_wtonlyoverlap <- findOverlaps(locAnn,methCHG_wtonly)
		methCHG_wtonlyoverlapunique <- !rev(duplicated(rev(queryHits(methCHG_wtonlyoverlap))))
		locAnn$methCHG_wtonly <- rep(FALSE, length(locAnn))
		locAnn$methCHG_wtonly[queryHits(methCHG_wtonlyoverlap)[methCHG_wtonlyoverlapunique]] <- TRUE
		locAnn$methCHG_wtonly <- ordered(locAnn$methCHG_wtonly)
 
		methCHG_dcl3onlyoverlap <- findOverlaps(locAnn,methCHG_dcl3only)
		methCHG_dcl3onlyoverlapunique <- !rev(duplicated(rev(queryHits(methCHG_dcl3onlyoverlap))))
		locAnn$methCHG_dcl3only <- rep(FALSE, length(locAnn))
		locAnn$methCHG_dcl3only[queryHits(methCHG_dcl3onlyoverlap)[methCHG_dcl3onlyoverlapunique]] <- TRUE
		locAnn$methCHG_dcl3only <- ordered(locAnn$methCHG_dcl3only)
        locAnn
}


# counting biases - uses the alignment data object used to run the segmentation
countingBiases <- function(locAnn, cl) {
    load("/home/bioinf/nem34/segmentation_with_externals.r_2015-11-12_17:56:17.079066/aDlt20_first_chlamy_segmentation_nick.RData") #aD #alignmentData
    colnames(values(aD@alignments))[2] <- "multireads"

    #samples <- read.delim("Arabidopsis_samples.txt", as.is = TRUE)
    
    samples <- read.csv("/home/bioinf/nem34/segmentMap_II/Summary_of_Data_5.csv")
    
    wt <- samples$GenuineControls %in% c("wt")
    
    aDwidths <- width(aD@alignments)
    aDnormal <- aD[aDwidths>19 & aDwidths < 22,] #Changed to more relevent for chlamy - was previously >21 and <24
    aDnormalnomulti <- aDnormal[!duplicated(as.character(aDnormal@alignments$tag)),]
    
#    save(aDnormal,file="aDnormal.rdata")
#    load("aDnormal.rdata")
    aDwidths <- width(aDnormal@alignments)
    
    firstNuc <- substr(aDnormal@alignments$tag,1,1)
    firstNucnomulti <- substr(aDnormalnomulti@alignments$tag,1,1)
                                        #proportion of sRNAs with a given 5'nuc
    table(firstNucnomulti)/length(firstNucnomulti)


    # find normal ratio of first base nucleotides
    expectedRatio <-  (tapply(rowSums(aDnormalnomulti@data),firstNucnomulti,sum)/sum(aDnormalnomulti@data))


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

    # assumes binomial distribution and looks for significant variation for each locus - variation from the overall ratio!
    z <- pbinom(firstBase,
                matrix(rowSums(firstBase), ncol = ncol(firstBase), nrow = nrow(firstBase)),
                prob = matrix(expectedRatio, ncol = ncol(firstBase), nrow = nrow(firstBase), byrow = TRUE), lower.tail = FALSE)
    z[rowSums(firstBase) == 0,] <- NA
    zq <- matrix(p.adjust(z, method = "BH"), ncol = ncol(firstBase))
    zq[rowSums(firstBase) == 0,] <- 1
    pred <- apply(do.call("cbind", lapply(1:4, function(ii) c("", c("A", "C", "G", "T")[ii])[as.integer(zq[,ii] < 0.01) + 1])), 1, paste, collapse = "")
    pred[pred == ""] <- NA
    
    locAnn$predominant_5prime_letter <- as.factor(pred)
    

                                        # counting numbers of 20s and 21s
    aDs <- aD[width(aD@alignments)<20,]
    countss <- getCounts(segments=locAnn,aD=aDs,cl=cl); rm(aDs)
    countsswt <- countss[,wt]
    aD20 <- aD[width(aD@alignments)==20,]
    counts20 <- getCounts(segments=locAnn,aD=aD20,cl=cl); rm(aD20)
    counts20wt <- counts20[,wt]
    aD21 <- aD[width(aD@alignments)==21]
    counts21 <- getCounts(segments=locAnn,aD=aD21,cl=cl); rm(aD21)
    counts21wt <- counts21[,wt]
	    aD21 <- aD[width(aD@alignments)==21]
    counts21 <- getCounts(segments=locAnn,aD=aD21,cl=cl); rm(aD21)
    counts21wt <- counts21[,wt]
	aDnorm <- aD[width(aD@alignments)<22 & width(aD@alignments)>19]
    countsnorm <- getCounts(segments=locAnn,aD=aDnorm,cl=cl); rm(aDnorm)
    countsnormwt <- countsnorm[,wt]
    aDb <- aD[width(aD@alignments)>21,]
    countsb <- getCounts(segments=locAnn,aD=aDb,cl=cl); rm(aDb)
    countsbwt <- countsb[,wt]
	
										#Counts for other sizes - Didn't yield anything interesting
	#aD22 <- aD[width(aD@alignments)==22]
    #counts22 <- getCounts(segments=locAnn,aD=aD22,cl=cl); rm(aD22)
    #counts22wt <- counts22[,wt]
	#aD23 <- aD[width(aD@alignments)==23]
    #counts23 <- getCounts(segments=locAnn,aD=aD23,cl=cl); rm(aD23)
    #counts23wt <- counts23[,wt]
	#aD24 <- aD[width(aD@alignments)==24]
    #counts24 <- getCounts(segments=locAnn,aD=aD24,cl=cl); rm(aD24)
    #counts24wt <- counts24[,wt]
	#aD25 <- aD[width(aD@alignments)==25]
    #counts25 <- getCounts(segments=locAnn,aD=aD25,cl=cl); rm(aD25)
    #counts25wt <- counts25[,wt]
	#aD26 <- aD[width(aD@alignments)==26]
    #counts26 <- getCounts(segments=locAnn,aD=aD26,cl=cl); rm(aD26)
    #counts26wt <- counts26[,wt]
	#aD27 <- aD[width(aD@alignments)==27]
    #counts27 <- getCounts(segments=locAnn,aD=aD27,cl=cl); rm(aD27)
    #counts27wt <- counts27[,wt]
	#aD28 <- aD[width(aD@alignments)==28]
    #counts28 <- getCounts(segments=locAnn,aD=aD28,cl=cl); rm(aD28)
    #counts28wt <- counts28[,wt]
	
	#locAnn$counts22 <- rowSums(counts22wt)
	#locAnn$counts23 <- rowSums(counts23wt)
    #locAnn$counts24 <- rowSums(counts24wt)
	#locAnn$counts25 <- rowSums(counts25wt)
    #locAnn$counts26 <- rowSums(counts26wt)
	#locAnn$counts27 <- rowSums(counts27wt)
    #locAnn$counts28 <- rowSums(counts28wt)
	
	#locAnn$ratio22vsNormal <- log2((rowSums(counts22wt))/(rowSums(countsnormwt)))
	#locAnn$ratio23vsNormal <- log2((rowSums(counts23wt))/(rowSums(countsnormwt)))
	#locAnn$ratio24vsNormal <- log2((rowSums(counts24wt))/(rowSums(countsnormwt)))
	#locAnn$ratio25vsNormal <- log2((rowSums(counts25wt))/(rowSums(countsnormwt)))
	#locAnn$ratio26vsNormal <- log2((rowSums(counts26wt))/(rowSums(countsnormwt)))
	#locAnn$ratio27vsNormal <- log2((rowSums(counts27wt))/(rowSums(countsnormwt)))
	#locAnn$ratio28vsNormal <- log2((rowSums(counts28wt))/(rowSums(countsnormwt)))
	
	#locAnn$ratio22vsNormalClass <- classCI(locAnn$counts22, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="22vsNormal_0.3_0.7.pdf")
	#locAnn$ratio23vsNormalClass <- classCI(locAnn$counts23, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="23vsNormal_0.3_0.7.pdf")
	#locAnn$ratio24vsNormalClass <- classCI(locAnn$counts24, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="24vsNormal_0.3_0.7.pdf")
	#locAnn$ratio25vsNormalClass <- classCI(locAnn$counts25, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="25vsNormal_0.3_0.7.pdf")
	#locAnn$ratio26vsNormalClass <- classCI(locAnn$counts26, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="26vsNormal_0.3_0.7.pdf")
	#locAnn$ratio27vsNormalClass <- classCI(locAnn$counts27, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="27vsNormal_0.3_0.7.pdf")
	#locAnn$ratio28vsNormalClass <- classCI(locAnn$counts28, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="28vsNormal_0.3_0.7.pdf")
        
                                        #sRNA size ratio
     # this function computes confidence intervals on the ratio of 20s and 21s for each locus, and then uses that to put each locus in a window. This function will produce a plot of the density of the mean (and the log of the mean) of the confidence intervals; choose the 'probs' values in such a way to split the modes of the density plots.

					#For 21 vs 20 ratio
	locAnn$ratio21vs20 <- log2((rowSums(counts21wt))/(rowSums(counts20wt)))
	locAnn$counts20 <- rowSums(counts20wt)
    locAnn$counts21 <- rowSums(counts21wt)
    locAnn$ratio21vs20[(locAnn$counts20+locAnn$counts21)<=5] <- NaN

    locAnn$ratio21vs20Class <- classCI(locAnn$counts21, locAnn$counts20, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="21vs20_0.3_0.7.pdf")

					#For 20 vs 21 ratio
    locAnn$ratio20vs21 <- log2((rowSums(counts20wt))/(rowSums(counts21wt)))

    locAnn$ratio20vs21[(locAnn$counts20+locAnn$counts21)<=5] <- NaN


    #locAnn$ratio20vs21Class <- classCI(locAnn$counts20, locAnn$counts21, probs = c(med = 10^-1.5, high = 10^-0.5), comma = 2,plotname="20vs21.pdf")
    #With different divisions
    locAnn$ratio20vs21Class <- classCI(locAnn$counts20, locAnn$counts21, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="20vs21_0.3_0.7.pdf")	
	
			#For small RNAs vs Normal RNAs
	locAnn$ratioSmallvsNormal <- log2((rowSums(countsswt))/(rowSums(countsnormwt)))
	locAnn$countsSmall <- rowSums(countsswt)
	locAnn$countsBig <- rowSums(countsbwt)
    locAnn$countsNormal <- rowSums(countsnormwt)
	locAnn$ratioSmallvsNormal <- log2((rowSums(countsswt))/(rowSums(countsnormwt)))
    locAnn$ratioSmallvsNormal[(locAnn$countsSmall+locAnn$countsNormal)<=5] <- NaN

    locAnn$ratioSmallvsNormalClass <- classCI(locAnn$countsSmall, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="SmallvsNormal_0.3_0.7.pdf")

				#For big RNAs vs Normal RNAs
	locAnn$ratioBigvsNormal <- log2((rowSums(countsbwt))/(rowSums(countsnormwt)))
    locAnn$ratioBigvsNormal[(locAnn$countsBig+locAnn$countsNormal)<=5] <- NaN
    locAnn$ratioBigvsNormalClass <- classCI(locAnn$countsBig, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="BigvsNormal_0.3_0.7.pdf")

				#For small RNAs vs Normal RNAs
	locAnn$ratioSmallvsNormal <- log2((rowSums(countsswt))/(rowSums(countsnormwt)))
    locAnn$ratioSmallvsNormal[(locAnn$countsSmall+locAnn$countsNormal)<=5] <- NaN
    locAnn$ratioSmallvsNormalClass <- classCI(locAnn$countsSmall, locAnn$countsNormal, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="SmallvsNormal_0.3_0.7.pdf")

				#For big RNAs vs small RNAs
	#locAnn$ratioBigvsSmall <- log2((rowSums(countsbwt))/(rowSums(countsswt)))
    #locAnn$ratioBigvsSmall[(locAnn$countsBig+locAnn$countsSmall)<=5] <- NaN
    #locAnn$ratioBigvsSmallClass <- classCI(locAnn$countsBig, locAnn$countsSmall, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="BigvsSmall_0.3_0.7.pdf")
	
				#For small RNAs bs Big RNAs
	#locAnn$ratioSmallvsBig <- log2((rowSums(countsswt))/(rowSums(countsbwt)))
    #locAnn$ratioSmallvsBig[(locAnn$countsSmall+locAnn$countsBig)<=5] <- NaN
    #locAnn$ratioSmallvsBigClass <- classCI(locAnn$countsSmall, locAnn$countsBig, probs = c(med = 0.3, high = 0.7), comma = 2,plotname="SmallvsBig_0.3_0.7.pdf")

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
                                         probs = c(0.1,0.35,0.65,0.9), comma = 1,plotname="plusvsminus2.pdf") #changed probs to 0.1 and 0.9 rather than 0.2 and 0.8
    levels(locAnn$ratio_strand_class) <- c("strong bias", "med bias", "no bias", "med bias", "strong bias")
    

        ##computing repetetivness for each loci (total reads div by multi match corrected) - use full setof aD not aDnormal
    matches <- aD@alignments$multireads
    countswt <- getCounts(segments=locAnn,aD=aD,cl=cl)
    aD@data <-  aD@data/matches
    
    countswtnorm <- getCounts(segments=locAnn,aD=aD,cl=cl)

    locAnn$repetitiveness <- 1-(rowSums(countswtnorm)/rowSums(countswt))
    locAnn$repetitivenessClass <- as.ordered(cut(locAnn$repetitiveness,unique(quantile(locAnn$repetitiveness,probs=seq(0,1,0.25),na.rm=TRUE),include.lowest=TRUE))) #Added unique function
	levels(locAnn$repetitivenessClass) <- c("low","median","high","very_high")
    locAnn
}

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

#Needed to compute confidence intervals in countingBiases function
classCI <- function(x1, x2, probs, divisions, comma, plot = TRUE,plotname) {
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
pdf(plotname)
    plot(density(log10(colMeans(ratioCI[,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))

    plot(density(colMeans(ratioCI[,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))
    p2r <- function(p) 1 / ((1 / p) - 1)
    plot(density(log2(p2r(colMeans(ratioCI[,!is.na(classIDs)]))))); abline(v = log2(p2r(probs)))
    plot(density(log10((ratioCI[1,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))
    plot(density((ratioCI[1,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))
    plot(density(log10((ratioCI[2,!is.na(classIDs)])), na.rm = TRUE)); abline(v = log10(probs))
    plot(density((ratioCI[2,!is.na(classIDs)]), na.rm = TRUE, from = 0, to = 1)); abline(v = (probs))

dev.off()
    
    classIDs
}

##Compares loci from different life cycles
lifeCycle <- function(locAnn, loci) {
    source("/home/tjh48/Code/segmentSeq_devel/segmentSeq/R/selectLoci.R")
    selLoc <- selectLoci(loci, FDR = 0.1, perReplicate = TRUE, returnBool = TRUE)
    samples <- read.csv("Summary_of_Data_5.csv",header=TRUE)
    #wt <- samples$Ecotype %in% c("wt")
#Control for genotype    
    CC1883 <- samples$Genotype %in% c("CC1883")

 
    locAnn$vegetative <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$LifeCycle == "Vegetative" & CC1883)]))]) > 0
    locAnn$zygote <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$LifeCycle == "Zygote" & CC1883)]))]) > 0
      
                                        
    locAnn$vegetativespecific <- locAnn$vegetative & !locAnn$zygote

    locAnn$zygotespecific <- locAnn$zygote & !locAnn$vegetative
    
    locAnn
}


##Compares loci from different strains
strainSpec <- function(locAnn, loci) {
    source("/home/tjh48/Code/segmentSeq_devel/segmentSeq/R/selectLoci.R")
    selLoc <- selectLoci(loci, FDR = 0.1, perReplicate = TRUE, returnBool = TRUE)
    samples <- read.csv("Summary_of_Data_5.csv",header=TRUE)
    wt <- samples$Ecotype %in% c("wt")

 
    locAnn$CC1883 <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Genotype == "CC1883" & wt)]))]) > 0
    locAnn$CC125 <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Genotype == "CC125" & wt)])),drop=FALSE]) > 0
    locAnn$J <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$Genotype == "J" & wt)])),drop=FALSE]) > 0
      
                                        
    locAnn$CC1883specific <- locAnn$CC1883 & !(locAnn$CC125 | locAnn$J)

    locAnn$CC125specific <- locAnn$CC125 & !(locAnn$CC1883 | locAnn$J)
    
    locAnn$Jspecific <- locAnn$J & !(locAnn$CC1883 | locAnn$CC125)
    locAnn
}

##Compares loci from mutant experiments - NOTE: selects only the wts from Adrian's mutant experiments for comparison!
mutantSpec <- function(locAnn, loci) {
    source("/home/tjh48/Code/segmentSeq_devel/segmentSeq/R/selectLoci.R")
    selLoc <- selectLoci(loci, FDR = 0.1, perReplicate = TRUE, returnBool = TRUE)
    samples <- read.csv("Summary_of_Data_5.csv",header=TRUE)
    #wt <- samples$Ecotype %in% c("wt")
#Control for genotype    
    #CC1883 <- samples$Genotype %in% c("CC1883")
 
    locAnn$wtAdrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$AdrianExp == "wt")]))]) > 0
    locAnn$dcl3Adrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(samples$AdrianExp == "dcl3")]))]) > 0
    locAnn$mutantAdrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(!samples$AdrianExp == "wt" & !samples$AdrianExp == "FALSE")]))]) > 0
    #locAnn$notdcl3Adrian <- rowSums(selLoc[,as.integer(unique(loci@replicates[which(!samples$AdrianExp == "FALSE" & !samples$AdrianExp == "dcl3")]))]) > 0

                                        #Specific to wild typeCC1883 from Adrian's experiments
    locAnn$wtAdrianspecific <- locAnn$wtAdrian & !locAnn$mutantAdrian
					#Specific to Adrian's mutants
    locAnn$mutantAdrianspecific <- locAnn$mutantAdrian & !locAnn$wtAdrian
    					
					#Specific to dcl3 related mutants
    locAnn$notDCL3dependent <- locAnn$dcl3Adrian & !locAnn$wtAdrian
					#Specifically not present in dcl3 mutants
    locAnn$DCL3dependent <- locAnn$wtAdrian & !locAnn$dcl3Adrian
    locAnn
}

#Including phasing outputs
phaseMatch <- function(locAnn)
    {
        ##TASI analyse

        beds <- as.data.frame(locAnn)[,1:3]
        #beds[,1] <- paste("Chr", beds[,1], sep = "")
        rownames(beds) <- 1:nrow(beds)
        write.table(beds, col.names = TRUE, row.names = TRUE, file = "/home/bioinf/nem34/segmentMap_II/src/phasing_loci.txt", quote = FALSE, sep = "\t")
        system("/home/bioinf/nem34/segmentMap_II/src/run_phasing2.py")
        
        tasi <- read.csv("/home/bioinf/nem34/segmentMap_II/phasing_results_by_locus_21nt.tsv",sep="\t",header=FALSE) 
                                        #table(locAnn$cluster[tasi[,"V10"]< c(-5)] & tasi[,"V4"]< c(-5)])
        Phased <- rep("none",nrow(tasi))
 
#If in more than 5 libraries loci phasing significant at 0.05 significance level = "high", between 1 and 5 = "moderate"
	tasi[,4:166] <- exp(tasi[,4:166])
        tasi[,4:166] <- apply(tasi[,4:166], 2, function(x) {x[x<1] <- p.adjust(x[x < 1], method = "BH"); x})
        Phased[rowSums(tasi[,4:166] < 0.05)>1] <- "moderate"
        Phased[rowSums(tasi[,4:166] < 0.05)>5] <- "high"


        
        tasi[,1] <- gsub("Chr", "", tasi[,1])
        tasmat <- match(apply(as.data.frame(locAnn)[,1:3], 1, function(x) paste(gsub(" ", "", x), collapse = ":")),
                        apply(as.data.frame(tasi)[,1:3], 1, function(x) paste(gsub(" ", "", x), collapse = ":")))
        locAnn$Phased <- "none"
        locAnn$Phased[!is.na(tasmat)] <- Phased[tasmat[!is.na(tasmat)]]        
        locAnn$Phased <- ordered(as.factor(locAnn$Phased), levels = c("none", "moderate", "high"))
        
        locAnn
    }





