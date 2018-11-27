
library(segmentSeq)
library(readr)
library(dplyr)
library(stringr)
baseDir      <- "/projects/nick_matthews"
segLocation  <- file.path(baseDir, "segmentation_2018")
gitdir       <- file.path(baseDir, "chlamy_locus_map_github")

# selected libraries as discussed with Nick/Seb/Adrian (july2018)
files <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
    filter(LociRun2018 == "Yes") %>%
    filter(Controls == "wt")

mylibs <- str_replace(basename(files$File)[1:2], "gff2.gz", "fasta")
libDir <- file.path(baseDir, "sequencing_data/fasta_all")
file.path(libDir, mylibs)
mylibs
# [1] "SL2108_L1.RUN540.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.fasta"
# [2] "SL2121_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.fast

tankcmd <- "perl PhaseTank_v1.0_mod.pl --genome /projects/nick_matthews/resources/Creinhardtii_236.fa --lib"

## needs finishing!!
