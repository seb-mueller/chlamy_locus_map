### Mass segementation of all available Chlamy data
#Author: nem34 (Using code from bacms2 previous segmentation)
#Date: 12/11/15


#Import Segment Seq
library(segmentSeq)

#Define data directory
datadir<-""


#Make cluster
cl<-makeCluster(46)

#define Chromosome lengths
chrlens <- c(8033585,9223677,9219486,4091191,3500558,9023763,6421821,5033832,7956127,6576019,
		3826814,9730733,5206065,4157777,1922860,7783580,7188315,271631,219038,200793,
		189560,163774,127913,127161,102191,80213,55320,55278,52813,52376,
		48183,42264,39192,33576,32450,25399,24537,24437,22408,22082,
		21325,21000,20974,17736,16939,16627,14746,14165,13462,
		12727,11225,6241,2479,2277)

#Define chromosomes and scaffolds
chrs <- c("chromosome_1","chromosome_2","chromosome_3","chromosome_4","chromosome_5","chromosome_6","chromosome_7","chromosome_8","chromosome_9","chromosome_10",
		"chromosome_11","chromosome_12","chromosome_13","chromosome_14","chromosome_15","chromosome_16","chromosome_17","scaffold_18","scaffold_19","scaffold_20",
		"scaffold_21","scaffold_22","scaffold_23","scaffold_24","scaffold_25","scaffold_26","scaffold_27","scaffold_28","scaffold_29","scaffold_30",
		"scaffold_31","scaffold_32","scaffold_33","scaffold_34","scaffold_35","scaffold_36","scaffold_37","scaffold_38","scaffold_39","scaffold_40",
		"scaffold_41","scaffold_42","scaffold_43","scaffold_44","scaffold_45","scaffold_46","scaffold_47","scaffold_48","scaffold_49","scaffold_50",
		"scaffold_51","scaffold_52","scaffold_53","scaffold_54")

#Define Columns in Input Data
cols = c(chr = 1, tag = 10, start = 4, end = 5, count = 6, strand = 7)

#Get Files
libfiles <- c(
	#Attila
		"/data/pipeline/prod/SL16/SL16.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
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
	#Andy - lots of 23/24
		"/data/pipeline/prod/SL165_1/SL165_1.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL165_2/SL165_2.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL166/SL166.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
	#Attila
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
	#Daisy
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
	#Attila
		"/data/pipeline/prod/SL14/SL14.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
	#Adrian
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
	#Betty
		"/data/pipeline/prod/SL2404/SL2404_L1.RUN626.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2405/SL2405_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2406/SL2406_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2407/SL2407_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2408/SL2408_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
		"/data/pipeline/prod/SL2409/SL2409_L1.RUNunknown.trimmed.filtered_15.non_redundant.v_genome_JGI_assembly5_Chlamydomonas_reinhardtii.patman.gff2",
	#External
		"/home/nem34/data_for_seg/GSM803103_CRE1_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/GSM803104_CRE2_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/GSM803105_CRE3_trimmed_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/Qi_Wang_sRNAs_aligned_assembly5.gff2",
		"/home/nem34/data_for_seg/schwach_sRNAs_aligned_assembly5.gff2")

#Name Libraries
libnames <- c(		
	#Attila	
		"SL16","SL17",	
		"SL60","SL74","SL75",	
		"SL76","SL139","SL140",	
		"SL141","SL142","SL143",	
		"SL144","SL145","SL146",
	#Andy		
		"SL165_1","SL165_2","SL166",
	#Attila
		"SL184","SL185","SL186",
		"SL187","SL254","SL255",
		"SL256","SL257","SL265",
		"SL304","SL305","SL306",
		"SL307","SL308",
		"SL379","SL380","SL390",
		"SL391","SL392","SL393",
	#Daisy
		"SL2108","SL2121","SL2122",
		"SL2124","SL2123","SL2125",
		"SL2126","SL2127","SL2128",
		"SL2129","SL2130","SL2131",
		"SL2132","SL2132","SL2134",
		"SL2135","SL2136","SL2137",
	#Attila
		"SL14",
	#Adrian
		"SL2178","SL2180","SL2179",
		"SL2181","SL2182","SL2183",
		"SL2184","SL2185","SL2186",
		"SL2187","SL2188","SL2189",
		"SL2190","SL2191","SL2192",
		"SL2193","SL2194","SL2195",
		"SL2196","SL2197","SL2198",
		"SL2199","SL2200","SL2201",
		"SL2202","SL2203","SL2204",
		"SL2205","SL2206","SL2207",
		"SL2208","SL2209","SL2210",
		"SL2211","SL2212","SL2213",
		"SL2214","SL2215","SL2216",
		"SL2217","SL2218","SL2219",
		"SL2220","SL2221","SL2222",
		"SL2223","SL2224","SL2225",
		"SL2298","SL2299","SL2300",
		"SL2301","SL2302","SL2303",
		"SL2304","SL2305","SL2306",
		"SL2307","SL2308","SL2309",
		"SL2310","SL2311","SL2312",
		"SL2313","SL2314","SL2315",
		"SL2316","SL2317","SL2318",
		"SL2319","SL2320","SL2321",
		"SL2322","SL2323","SL2324",
		"SL2325","SL2326","SL2327",
		"SL2328","SL2329","SL2330",
		"SL2331","SL2332","SL2333",
		"SL2362","SL2363","SL2364",
		"SL2365","SL2366","SL2367",
		"SL2368","SL2369","SL2370",
		"SL2371","SL2372","SL2373",
	#Betty	
		"SL2404","SL2405","SL2406",
		"SL2407","SL2408","SL2409",
	#External
		"GSM803103","GSM803104","GSM803105",
		"GSM183546","GSM176482")


#Define replicates
replicates<-as.factor(c(
		#Attila
			1,2,
			3,4,5,
			6,7,7,
			8,8,9,
			9,10,10,
		#Andy
			11,11,12,
		#Attila
			13,14,15,
			16,17,18,
			19,20,21,
			22,23,24,
			25,26,
			27,28,29,
			29,30,30,
		#Daisy
			31,31,31,
			32,32,32,
			33,33,33,
			34,34,34,
			35,35,35,
			36,36,36,
		#Attila
			37,
		#Adrian
			38,38,38,
			39,39,39,
			40,40,40,
			41,41,41,
			42,42,42,
			43,43,43,
			44,44,44,
			45,45,45,
			46,46,46,
			47,47,47,
			48,48,48,
			49,49,49,
			50,50,50,
			51,51,51,
			52,52,52,
			53,53,53,
			54,54,54,
			55,55,55,
			56,56,56,
			57,57,57,
			58,58,58,
			59,59,59,
			60,60,60,
			61,61,61,
			62,62,62,
			63,63,63,
			64,64,64,
			65,65,65,
			66,66,66,
			67,67,67,
			68,68,68,
			69,69,69,
		#Betty
			70,70,70,
			71,71,71,
		#External
			72,73,74,
			75,76))



#creating alignment objext
aD <- readGeneric(files = libfiles,dir = datadir,replicates = replicates, libnames = libnames,chrs = chrs, chrlens = chrlens,cols=cols, verbose=TRUE, gap = 200,cl=cl)
save(aD, file="aD_first_chlamy_segmentation_nick.RData")
#load("aD_first_chlamy_segmentation_nick.RData")

#Get rid of highly expressed and repetative data
aD<-aD[aD@alignments$multireads<20]
save(aD, file="aDlt20_first_chlamy_segmentation_nick.RData")
#load("aDlt20_first_chlamy_segmentation_nick.RData")

#Process alignment data to find potential segements
sD<-processAD(aD,gap=200, cl=cl) #How big a gap should I use?
save(sD, file="sDlt20_first_chlamy_segmentation_nick.RData")
#load("sDlt20_first_chlamy_segmentation_nick.RData")

#Generate Locus map
hS<-heuristicSeg(sD=sD,aD=aD,getLikes=TRUE,cl=cl)
save(hS, file="hS_first_chlamy_segmentation_nick.RData")
#load("hS_first_chlamy_segmentation_nick.RData")0

#Generate a genome map
segD<-classifySeg(aD=aD,sD=sD,cD=hS, getLikes=TRUE,cl=cl)
save(segD, file="segD_first_chlamy_segmentation_nick.RData")
#load("segD_first_chlamy_segmentation_nick.RData")
