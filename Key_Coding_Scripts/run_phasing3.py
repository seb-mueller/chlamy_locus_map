#Script written by Bruno Santos to analyse phasing of loci
#Date: 01/12/15

#!/usr/bin/env python
"""
24 Sep 2012
Python Script to perform search of Phased sRNAs
on the libraries used to test the algorithm
It takes as input a list of loci, collects the
sRNAs mapping to them and finally runs the phasing
method on the locus
 """
from Methods import Methods 

#Standard call to improve speed 
methods = Methods()
phaser = methods.phaser
chen = methods.chen
brachy = methods.brachy


class RNAsignature():
    def __init__(self,zChr,nCoord,nStrand,zTag,nCounts,nsize,nEdits=0):
        """
        RNAsignature(self,zChr,nCoord,nStrand,zTag,nCounts,nsize)
        Store the information about each read in an alignment file
        zChr <- Target of the alignment
        nCoord <- Start position of the read
        nStrand <- Direct or inverted read
        zTag <- Nucleotide Sequence
        nCounts <- Number of tags from the sequencing
        """
        self.zChr = zChr
        self.nCoord = int(nCoord)
        self.nStrand = int(nStrand)
        self.zTag = zTag
        self.nCounts = float(nCounts)
        self.nsize = int(nsize)
        self.nEdits = nEdits
        
if __name__=='__main__':
    #Define location of the pipeline processed files
    directory = '' #Directory containing the gff2 files
    lslibs = ["/data/pipeline/prod/SL16/SL16.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL17/SL17.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL60/SL60.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL74/SL74.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",	
		"/data/pipeline/prod/SL75/SL75.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL76/SL76.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL139/SL139.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL140/SL140.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL141/SL141.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL142/SL142.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL143/SL143.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL144/SL144.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL145/SL145.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL146/SL146.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL165_1/SL165_1.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL165_2/SL165_2.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL166/SL166.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL184/SL184.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL185/SL185.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL186/SL186.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL187/SL187.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL254/SL254.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL255/SL255.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL256/SL256.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL257/SL257.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL265/SL265.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL304/SL304.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL305/SL305.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL306/SL306.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL307/SL307.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL308/SL308.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL379/SL379.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL380/SL380.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL390/SL390.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL391/SL391.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL392/SL392.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL393/SL393.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2108/SL2108_L1.RUN540.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2121/SL2121_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2122/SL2122_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2124/SL2124_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2123/SL2123_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2125/SL2125_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2126/SL2126_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2127/SL2127_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2128/SL2128_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2129/SL2129_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2130/SL2130_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2131/SL2131_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2132/SL2132_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2133/SL2133_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2134/SL2134_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2135/SL2135_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2136/SL2136_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2137/SL2137_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL14/SL14.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2178/SL2178_L1.RUN546.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2180/SL2180_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",		
		"/data/pipeline/prod/SL2179/SL2179_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2181/SL2181_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2182/SL2182_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2183/SL2183_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2184/SL2184_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2185/SL2185_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2186/SL2186_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2187/SL2187_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2188/SL2188_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2189/SL2189_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2190/SL2190_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2191/SL2191_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2192/SL2192_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2193/SL2193_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2194/SL2194_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2195/SL2195_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2196/SL2196_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2197/SL2197_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2198/SL2198_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2199/SL2199_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2200/SL2200_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2201/SL2201_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2202/SL2202_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2203/SL2203_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2204/SL2204_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2205/SL2205_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2206/SL2206_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2207/SL2207_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2208/SL2208_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2209/SL2209_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2210/SL2210_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2211/SL2211_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2212/SL2212_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2213/SL2213_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2214/SL2214_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2215/SL2215_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2216/SL2216_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2217/SL2217_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2218/SL2218_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2219/SL2219_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2220/SL2220_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2221/SL2221_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2222/SL2222_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2223/SL2223_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2224/SL2224_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2225/SL2225_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2298/SL2298_L1.RUN618.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2299/SL2299_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2300/SL2300_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2301/SL2301_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2302/SL2302_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2303/SL2303_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2304/SL2304_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2305/SL2305_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2306/SL2306_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2307/SL2307_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2308/SL2308_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2309/SL2309_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2310/SL2310_L1.RUN619.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2311/SL2311_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2312/SL2312_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2313/SL2313_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2314/SL2314_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2315/SL2315_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2316/SL2316_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2317/SL2317_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2318/SL2318_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2319/SL2319_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2320/SL2320_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2321/SL2321_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2322/SL2322_L1.RUN620.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2323/SL2323_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2324/SL2324_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2325/SL2325_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2326/SL2326_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2327/SL2327_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2328/SL2328_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2329/SL2329_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2330/SL2330_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2331/SL2331_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2332/SL2332_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2333/SL2333_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2362/SL2362_L1.RUN622.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2363/SL2363_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2364/SL2364_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2365/SL2365_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2366/SL2366_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2367/SL2367_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2368/SL2368_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2369/SL2369_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2370/SL2370_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2371/SL2371_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2372/SL2372_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2373/SL2373_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2404/SL2404_L1.RUN626.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2405/SL2405_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2406/SL2406_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2407/SL2407_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2408/SL2408_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2409/SL2409_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/home/nem34/data_for_seg/GSM803103_CRE1_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/GSM803104_CRE2_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/GSM803105_CRE3_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/Qi_Wang_sRNAs_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/schwach_sRNAs_aligned_assembly5.gff2"]
    
    #Define the size of sRNAs
    nsize = 21
    
    #Define location of the file containing the loci definition
    zname = '/home/bioinf/nem34/segmentMap_II/lociforphasing.csv'
    loci_handle = open(zname)
    print 'loci file being used is:%s' % zname
    #Read the loci file from segmentSeq
    dtloci = {}
    first_line = True
    for zline in loci_handle:
        if first_line:
            first_line=False
            continue
        nlnumber,zchr,zstart,zend,zwidth,zstrand = zline.rstrip('\n').split(',')      
        try:dtloci[zchr.strip('"')].append((int(zstart),int(zend)))
        except KeyError: dtloci[zchr.strip('"')] = [(int(zstart),int(zend))]
    loci_handle.close()
    
    print 'Loci file successfully processed'

    
    #Sort the list of loci per start coordinate
    for zchr in dtloci.iterkeys():dtloci[zchr].sort(key=lambda nstart:nstart[0])
    
    #Read in the gff files
    #Store the srnas in each library by chr and start coordinate
    print 'Reading in gff files'
    ls_srnas = []
    for zfile in lslibs:
        print zfile
        zname = zfile.split('/')[0]
        handle = open(directory+zfile)
        dtsignatures = {}
        dtConvert = {'+':'1','-':'-1'}
        for zline in handle:
            #1    SL258    patman    3    21    1    +    .    SL258 CTAAACCCTAAACCCTAAA
            zline = zline.rstrip('\n')
            zchr,zlib,zprogram,zstart,zend,zcount,zstrand,none,ztag = zline.split('\t')
            
            if int(zend)-int(zstart)+1!=nsize:continue
            if zstrand=='-':
                ztemp = zend
                zend = zstart
                zstart = ztemp
            try:dtsignatures[zchr].append((int(zstart),RNAsignature(zchr,zstart,dtConvert[zstrand],ztag.split()[-1],zcount,nsize)))    
            except KeyError: dtsignatures[zchr] = [(int(zstart),RNAsignature(zchr,zstart,dtConvert[zstrand],ztag.split()[-1],zcount,nsize))]
        #Close the file
        handle.close()

	#Sort dtsignatures by start coordinate of the sRNA
        for lssignatures in dtsignatures.itervalues():lssignatures.sort(key = lambda zstart:zstart[0])
        
        #Store the dtsignatures in the lslibs
        ls_srnas.append(dtsignatures)
        

    
    print 'library files successfully processed'
        
    #Find overlaps between loci and sRNAs
    
    #Record the dictionaries by locus
    dt_signatures_by_loci = {}
    #Iterate over each chromosome
    for zchr in dtloci.iterkeys():
        print 'Now processing chromosome %s' % zchr
        #Get the list of loci on that chromossome
        lsloci = dtloci[zchr]
        #Indicate the position in the list for each library
        pos_ls = [0]*len(lslibs)
        #Iterate over each locus
        loci_number = 0
        for nstart,nend in lsloci:
            #Print progress
            #if loci_number%100 == 0 :print 'now processing loci %d of %d' % (loci_number,len(lsloci))
            loci_number+=1
            #Start a new dtsignatures dictionary that will contain all the RNASignatures objects for that locus
            #dtsignatures = {}
            #Iterate over all the libraries
            lib_pos = 0 # Library location
            for dtsrnas in ls_srnas:
                #Start a new dtsignatures dictionary that will contain all the RNASignatures objects for that locus
                dtsignatures = {}
                try:lssignatures = dtsrnas[zchr]
                except KeyError:
                    try:dt_signatures_by_loci[zchr+','+str(nstart)+','+str(nend)].append(dtsignatures)
                    except KeyError: dt_signatures_by_loci[zchr+','+str(nstart)+','+str(nend)] = [dtsignatures]
                    continue
                sig_coord,signature = lssignatures[pos_ls[lib_pos]]
                while sig_coord <= nend and pos_ls[lib_pos]+1<len(lssignatures):
                    if sig_coord>=nstart:dtsignatures[signature.zChr+','+str(signature.nCoord)+','+str(signature.nStrand)]=signature
                    pos_ls[lib_pos]+=1
                    sig_coord,signature = lssignatures[pos_ls[lib_pos]]
                try:dt_signatures_by_loci[zchr+','+str(nstart)+','+str(nend)].append(dtsignatures)
                except KeyError: dt_signatures_by_loci[zchr+','+str(nstart)+','+str(nend)] = [dtsignatures]
                lib_pos+=1
        
    handle = open('/home/nem34/segmentMap_II/phasing_results_by_locus_%dnt.tsv' % nsize,'w')
    
    #Perform phasing analysis of the locus
    total = len(dt_signatures_by_loci)
    i = 0
    for zname,locus in dt_signatures_by_loci.iteritems():
        if i%100 == 0 :print 'now processing loci %d of %d' % (i,total)
        handle.write('\t'.join(zname.split(',')))
        for dtsignatures in locus:
            fpvalue,zchr,zstart,zend,nCounts = phaser(dtsignatures,nsize,0,10*nsize,10*nsize)#Parameters for phaser are srnas_dictionary,size,min_pValue,min_length,max_gap
            handle.write('\t%s' % (fpvalue))
            #fpvalue,zchr,zstart,zend =brachy(dtsignatures,nsize,10)#Parameters for Brachypodium method are: srnas_dictionary,size,number_of_cycles
            #handle.write('\t%s\t%s' % (fpvalue,'_'.join([zchr,zstart,zend])))
        handle.write('\n')
        i+=1
    handle.close()
            
