# chlamy_locus_map
Small RNA Locus Map for Chlamydomonas reinhardtii

## Paper location
https://www.overleaf.com/17574001mfpjmkhghhtc

## Important scripts/files in project:

### general_setup.sh

### chlamy_source_code.R
Source code file containing all functions applied in other scripts

### chlamy_segmentation_gap100.R
Mass segementation of all available Chlamy data
```
in:
"Summary_of_Data.csv" $File
out:
save(aD, meta, file=file.path(segLocation ,paste0("aD_chlamy_segmentation_",mycomment,".RData")))
save(hS, meta, file=file.path(segLocation, paste0("hS_chlamy_segmentation_", mycomment, ".RData")))
save(segD, meta, file=file.path(segLocation ,paste0("segD_chlamy_segmentation_", mycomment, ".RData")))
```

### Segmentation_Analysis.R
Analyses quality of segmentation producing a variety of graphs
```
in: "segD_chlamy_segmentation_multi200_gap100.RData"
out: only plots??
```

### Segmentation_Analysis_post_selection.R
 chlamy small RNA loci R code
 This code is based on a finished loci definition as obtained from:
 chlamy_annotation_pipeline.R and Segmentation_Analysis.R
 used to investigate loci for interesting association etc and descriptive features for writing paper
```
in:
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
```
### chlamy_annotation_pipeline.R
Script for running annotation functions as well as defining sRNA loci for chlamy.
 The appropriate thresholds (e.g. FDR) are determined in:
 Segmentation_Analysis.R
```
in:
gitdir, "Summary_of_Data.csv"
source(file.path(gitdir, "Scripts/chlamy_source_code.R"))
inputdata <- "segD_chlamy_segmentation_LociRun2018_multi200_gap100.RData" #14390
aDfile    <- "aD_chlamy_segmentation_LociRun2018_multi200_gap100.RData"
"transposon_annotations.Rdata"
out:
save(gr, meta, metawt, loci, baseDir, prefix, saveLocation, file = file.path(saveLocation, paste0("gr_fdr", fdr,  ".RData")))
export.gff3(gr, con = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".gff")))
write.csv(as.data.frame(gr), file = file.path(saveLocation, paste0("loci_fdr", fdr, prefix, ".csv")))
```

### chlamy_MCA_analysisis.R
work out which clusters and dimensions to use, and takes ages
Script to compute computes diagnostic plots and figures for MCA and clustering
```
save(dimList, file = file.path(saveLocation,"dimList.RData"))
save(dimStab, file = file.path(saveLocation,"dimStab.RData"))
```
### chlamy_MCA_analysis_v2.R
 stability analysis etc. which takes a while (~2 hours)

### chlamy_MCA_analysis_plotting.R
takes the analysis output and plots them. 

###  chlamy_MCA.R
performs and plots the final MCA and HCPC based on the settings determined by the previous analysis. You'll find all the latest outputs and figures in "LociRun2018_multi200_gap100_90c7213_MCAOutputs_05c5bb5" including the heatmaps and chromosome tracks.
Script to run MCA to cluster loci according to their annotations
```
in:
gitdir,"Annotation2Use.csv"
inputFile <- "gr_fdr0.05.RData"
lociRun <- "LociRun2018_multi200_gap100_90c7213"
#baseDir <- "C:/Users/Nick/Documents/PhD/Projects/Chlamy"
# Specify location of annotation outputs
inputLocation <- file.path(baseDir, "segmentation_2018", MCAOutputs)
lociLocation <- file.path(baseDir, "segmentation_2018", lociRun)
out: resMCA, gr
MCAOutputs <- "LociRun2018_multi200_gap100_90c7213_MCAOutputs_05c5bb5"
```

### chlamy_summary_plots.R
Script to run MCA to cluster loci according to their annotations
Specify clusters and dimensions used
```
in:
load(file.path(annoDir,"chlamy_all_annotations.Rdata"), verbose = FALSE)
load(file.path(inputLocation,"clusterings.RData"))
load(file.path(inputLocation,"gr_clustered.RData"))
load(file.path(inputLocation,"resMCA.RData"))
nclust <- 6; ndim <- 7
out:
write.csv(chrLociTable,file=file.path(figLocation,"chrLociTable.csv"))
write.csv(annotTable,file=file.path(figLocation,"annotTable.csv"))
write.csv(clusterAnnotTable,file=file.path(figLocation,"clusterAnnotTable.csv"))
Clustercoverage.png
```
