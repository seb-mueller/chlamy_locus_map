
library(segmentSeq)
library(readr)
library(dplyr)
library(stringr)
baseDir      <- "/projects/nick_matthews"
segLocation  <- file.path(baseDir, "segmentation_2018")
gitdir       <- file.path(baseDir, "chlamy_locus_map_github")
setwd(file.path(gitdir, "PhaseTank_v1.0"))

# selected libraries as discussed with Nick/Seb/Adrian (july2018)
files <- read_csv(file.path(gitdir, "Summary_of_Data.csv")) %>%
    filter(LociRun2018 == "Yes") %>%
    filter(Controls == "wt")

files$DataCode
# [1] "SL2108" "SL2121" "SL2122" "SL2123" "SL2124" "SL2125"
# [7] "SL2181" "SL2182" "SL2183" "SL2184" "SL2185" "SL2186"
# [13] "SL2187" "SL2188" "SL2189" "SL2301" "SL2302" "SL2303"
# [19] "SL2310" "SL2311" "SL2312" "SL2313" "SL2314" "SL2315"
# [25] "SL2322" "SL2323" "SL2324" "SL2325" "SL2326" "SL2327"

libDir <- file.path(baseDir, "sequencing_data/fasta_with_counts")
mylibs <- file.path(libDir, paste0(files$DataCode,"_assembly5_Chlamydomonas_reinhardtii.patman.aligned_reads.fasta"))
head(mylibs)
# [1] "/projects/nick_matthews/sequencing_data/fasta_with_counts/SL2108_assembly5_Chlamydomonas_reinhardtii.patman.aligned_reads.  fasta"
# [2] "/projects/nick_matthews/sequencing_data/fasta_with_counts/SL2121_assembly5_Chlamydomonas_reinhardtii.patman.aligned_reads.  fasta"
file.exists(mylibs)

tankcmd <- "perl PhaseTank_v1.0_mod.pl --genome /projects/nick_matthews/resources/Creinhardtii_236.fa --lib"

system(paste(tankcmd, paste(mylibs, collapse = ",")))

#You can check the output files in directory './OUTPUT_2018.11.27_18.08/'.
