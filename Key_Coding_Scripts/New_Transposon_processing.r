#Code used to filter down transposons and have consistent naming scheme
#Repeat masked file used for annotation, code written by Nick Matthews
#Date: 21/01/15

library(rtracklayer)
annot<-import.gff3("/home/nem34/Creinhardtii_281_v5.5.repeatmasked_assembly_v5.0.gff3")

x<-grep("rich",annot$Name)
annot<-annot[-x]
x<-grep(")n",annot$Name)
annot<-annot[-x]
x<-grep("rnd",annot$Name)
annot<-annot[-x]

y<-grep("Gypsy",annot$Name)
 table(annot$Name[y])


unique(annot$Name)
  [1] "L1-6_CR"          "Gypsy8-I_CR-int"  "hAT-N2_CR"
  [4] "TE1-2_CR"         "hAT-N5_CR"        "Gypsy-5-I_CR-int"
  [7] "Gypsy11-I_CR-int" "DNA-1-7_CR"       "NonLTR-4_CR"
 [10] "TE2-7a_CR"        "REP2_CR"          "REP2A_CR"
 [13] "RTE-3_CR"         "NonLTR-1_CR"      "hAT-N11a_CR"
 [16] "Gypsy-6-LTR_CR"   "TE2auto_CR"       "L1-3_CR"
 [19] "TE2-5_CR"         "Copia2-I_CR-int"  "Novosib-N1"
 [22] "P-2_CR"           "Gypsy-3-LTR_CR"   "L1-2_CR"
 [25] "hAT-N11_CR"       "Novosib-3_CR"     "REPX2_CR"
 [28] "MarinerL-1_CR"    "hAT-N4_CR"        "hAT-N6_CR"
 [31] "Gypsy-4-I_CR-int" "L1-1_CR"          "Gypsy-4-LTR_CR"
 [34] "SINEX-4_CR"       "hAT-N9_CR"        "RandI-1"
 [37] "TE2-7b_CR"        "REPY1_CR"         "SINEX-3_CR"
 [40] "RandI-2"          "Mariner-N4_CR"    "DNA-4-7_CR"
 [43] "Gypsy-5-LTR_CR"   "LTR1_CR"          "RTE-2a_CR"
 [46] "TE1-1_CR"         "P-1_CR"           "DNA-8-1_CR"
 [49] "RTE-1_CR"         "Gypsy-1-I_CR-int" "Copia1-I_CR-int"
 [52] "hAT-N7_CR"        "DNAX1_CR"         "Copia1-LTR_CR"
 [55] "Gypsy9-LTR_CR"    "Gypsy9-I_CR-int"  "Gypsy11-LTR_CR"
 [58] "RTE-2_CR"         "DNA-3-7_CR"       "DIRS-1_CR"
 [61] "EnSpm-N1a_CR"     "DNA-2-7_CR"       "RandI-4"
 [64] "DNA3-1_CR"        "RandI-5"          "RandI-3"
 [67] "Novosib-2_CR"     "L1-4_CR"          "Mariner-N6_CR"
 [70] "L1-5_CR"          "TOC2"             "CRTOC1"
 [73] "Gulliver"         "EnSpm-N3_CR"      "nonLTR-5_CR"
 [76] "hAT-N3_CR"        "EnSpm-N2_CR"      "RandI-4A_CR"
 [79] "LTR2_CR"          "TE2-4_CR"         "TE2-1_CR"
 [82] "RTE-4_CR"         "TE2-2_CR"         "RandI-4B_CR"
 [85] "DNA-1-9_CR"       "MSAT-1E_CR"       "Novosib-N1a"
 [88] "Copia3-I_CR-int"  "TE2-3_CR"         "Copia3-LTR_CR"
 [91] "polypyrimidine"   "NonLTR-3_CR"      "Mariner-N5_CR"
 [94] "DNA-8-5_CR"       "DNA-8-4_CR"       "Gypsy-5B-LTR_CR"
 [97] "RandI-6"          "SUBTEL_sa"        "LTR3_CR"
[100] "TE2-6_CR"         "TCR1_LTR"         "polypurine"
[103] "REP1_CR"          "REP4_CR"          "EnSpm-N1_CR"
[106] "REM1-LTR"         "Gypsy-5A-LTR_CR"  "REPX4_CR"
[109] "Gypsy12-I_CR-int" "REPX3_CR"         "Gypsy-6-I_CR-int"
[112] "MSAT-1D_CR"       "hAT-N10_CR"       "Gypsy8-LTR_CR"
[115] "DNA-8-2_CR"       "hAT-N1_CR"        "SINEX-1_CR"
[118] "Gypsy-1-LTR_CR"   "Novosib-1_CR"     "hAT-N8_CR"
[121] "Gypsy-2-LTR_CR"   "RTE-5_CR"         "REM1-I-int"
[124] "Gypsy12-LTR_CR"   "SINEX-2_CR"       "MSAT-4A_CR"
[127] "Mariner-N1_CR"    "NonLTR-2_CR"      "Gypsy-2-I_CR-int"
[130] "MSAT-2_CR"        "Mariner-N3_CR"    "DNA-8-3_CR"
[133] "DNA-8-6_CR"       "REPX1_CR"         "Gypsy10-LTR_CR"
[136] "Gypsy10-I_CR-int" "Copia2-LTR_CR"    "Mariner-N2_CR"
[139] "Novosib-4_CR"

x<-grep("SUBTEL",annot$Name)
annot<-annot[-x]
x<-grep("MSAT",annot$Name)
annot<-annot[-x]
x<-grep("poly",annot$Name)
annot<-annot[-x]


export.gff3(annot,"/home/nem34/annot_newtrans.gff3")
#annot<-import.gff3("/home/nem34/annot_newtrans.gff3")

x<-grep("rnd",annot$Name)
annot<-annot[-x]

x<-grep("Gypsy",annot$Name)
annot$Name[x]<-"LTR.Gypsy"
x<-grep("Copia",annot$Name)
annot$Name[x]<-"LTR.Copia" 


  [1] "REP2_CR"        "LTR2_CR"        "NonLTR-1_CR"    "L1-6_CR"
  [5] "LTR.Gypsy"      "LTR.Copia"      "REP2A_CR"       "hAT-N2_CR"
  [9] "TE1-2_CR"       "hAT-N5_CR"      "RTE-2a_CR"      "RTE-4_CR"
 [13] "hAT-N11a_CR"    "DNA-1-7_CR"     "NonLTR-4_CR"    "hAT-N6_CR"
 [17] "TE2-7a_CR"      "DNA-2-7_CR"     "Novosib-3_CR"   "RTE-3_CR"
 [21] "hAT-N8_CR"      "L1-4_CR"        "L1-3_CR"        "TE2auto_CR"
 [25] "L1-5_CR"        "TE2-5_CR"       "Novosib-N1"     "RandI-4"
 [29] "P-2_CR"         "hAT-N10_CR"     "RTE-1_CR"       "L1-2_CR"
 [33] "hAT-N11_CR"     "hAT-N3_CR"      "REPX2_CR"       "hAT-N7_CR"
 [37] "MarinerL-1_CR"  "LTR4_CR"        "nonLTR-5_CR"    "hAT-N4_CR"
 [41] "hAT-N9_CR"      "REM1-I-int"     "RandI-6"        "RandI-1"
 [45] "L1-1_CR"        "SINEX-4_CR"     "TE2-7b_CR"      "REPY1_CR"
 [49] "SINEX-3_CR"     "RandI-2"        "RandI-3"        "EnSpm-N2_CR"
 [53] "Mariner-N4_CR"  "DNA-4-7_CR"     "REP4_CR"        "RTE-2_CR"
 [57] "LTR1_CR"        "TE2-3_CR"       "TE1-1_CR"       "P-1_CR"
 [61] "DNA-8-1_CR"     "DNA3-1_CR"      "DNAX1_CR"       "DNA-8-4_CR"
 [65] "DNA-3-7_CR"     "DIRS-1_CR"      "EnSpm-N1a_CR"   "SINEX-1_CR"
 [69] "RandI-5"        "Novosib-2_CR"   "Mariner-N6_CR"  "TOC2"
 [73] "CRTOC1"         "Gulliver"       "EnSpm-N3_CR"    "TE2-1_CR"
 [77] "RandI-4A_CR"    "TE2-2_CR"       "TE2-4_CR"       "DNA-1-9_CR"
 [81] "DNA-8-5_CR"     "DNA-8-2_CR"     "Novosib-N1a"    "REPX4_CR"
 [85] "Mariner-N2_CR"  "RandI-4B_CR"    "REP1_CR"        "hAT-N1_CR"
 [89] "Novosib-4_CR"   "NonLTR-3_CR"    "Mariner-N5_CR"  "NonLTR-2_CR"
 [93] "LTR3_CR"        "TE2-6_CR"       "TCR1_LTR"       "REM1-LTR"
 [97] "EnSpm-N1_CR"    "Novosib-1_CR"   "RTE-5_CR"       "REPX3_CR"
[101] "DNA-8-6_CR"     "SINEX-2_CR"     "Mariner-N3_CR"  "Mariner-N1_CR"
[105] "DNA-8-3_CR"     "REPX1_CR"       "Harbinger-1_CR"

x<-grep("Gypsy",annot$Name)
annot$Name[x]<-"Gypsy"
x<-grep("Copia",annot$Name)
annot$Name[x]<-"Copia" 
x<-grep("Novosib",annot$Name)
annot$Name[x]<-"Novosib" 
x<-grep("Mariner",annot$Name)
annot$Name[x]<-"Mariner" 
x<-grep("hAT",annot$Name)
annot$Name[x]<-"hAT" 
x<-grep("RTE",annot$Name)
annot$Name[x]<-"RTE" 
x<-grep("TOC1",annot$Name)
annot$Name[x]<-"TOC1" 
x<-grep("TOC2",annot$Name)
annot$Name[x]<-"TOC2" 

#10968 remaining
unique(annot$Name)

 [1] "REP2_CR"        "LTR2_CR"        "NonLTR-1_CR"    "L1-6_CR"
 [5] "Gypsy"          "Copia"          "REP2A_CR"       "hAT"
 [9] "TE1-2_CR"       "RTE"            "DNA-1-7_CR"     "NonLTR-4_CR"
[13] "TE2-7a_CR"      "DNA-2-7_CR"     "Novosib"        "L1-4_CR"
[17] "L1-3_CR"        "TE2auto_CR"     "L1-5_CR"        "TE2-5_CR"
[21] "RandI-4"        "P-2_CR"         "L1-2_CR"        "REPX2_CR"
[25] "Mariner"        "LTR4_CR"        "nonLTR-5_CR"    "REM1-I-int"
[29] "RandI-6"        "RandI-1"        "L1-1_CR"        "SINEX-4_CR"
[33] "TE2-7b_CR"      "REPY1_CR"       "SINEX-3_CR"     "RandI-2"
[37] "RandI-3"        "EnSpm-N2_CR"    "DNA-4-7_CR"     "REP4_CR"
[41] "LTR1_CR"        "TE2-3_CR"       "TE1-1_CR"       "P-1_CR"
[45] "DNA-8-1_CR"     "DNA3-1_CR"      "DNAX1_CR"       "DNA-8-4_CR"
[49] "DNA-3-7_CR"     "DIRS-1_CR"      "EnSpm-N1a_CR"   "SINEX-1_CR"
[53] "RandI-5"        "TOC2"           "TOC1"           "Gulliver"
[57] "EnSpm-N3_CR"    "TE2-1_CR"       "RandI-4A_CR"    "TE2-2_CR"
[61] "TE2-4_CR"       "DNA-1-9_CR"     "DNA-8-5_CR"     "DNA-8-2_CR"
[65] "REPX4_CR"       "RandI-4B_CR"    "REP1_CR"        "NonLTR-3_CR"
[69] "NonLTR-2_CR"    "LTR3_CR"        "TE2-6_CR"       "TCR1_LTR"
[73] "REM1-LTR"       "EnSpm-N1_CR"    "REPX3_CR"       "DNA-8-6_CR"
[77] "SINEX-2_CR"     "DNA-8-3_CR"     "REPX1_CR"       "Harbinger-1_CR"

x<-grep("DIRS",annot$Name)
annot$Name[x]<-"DIRS"
x<-grep("Harbinger",annot$Name)
annot$Name[x]<-"Harbinger" 
x<-grep("SINE",annot$Name)
annot$Name[x]<-"SINE"

unique(annot$Name)

 [1] "REP2_CR"      "LTR2_CR"      "NonLTR-1_CR"  "L1-6_CR"      "Gypsy"
 [6] "Copia"        "REP2A_CR"     "hAT"          "TE1-2_CR"     "RTE"
[11] "DNA-1-7_CR"   "NonLTR-4_CR"  "TE2-7a_CR"    "DNA-2-7_CR"   "Novosib"
[16] "L1-4_CR"      "L1-3_CR"      "TE2auto_CR"   "L1-5_CR"      "TE2-5_CR"
[21] "RandI-4"      "P-2_CR"       "L1-2_CR"      "REPX2_CR"     "Mariner"
[26] "LTR4_CR"      "nonLTR-5_CR"  "REM1-I-int"   "RandI-6"      "RandI-1"
[31] "L1-1_CR"      "SINE"         "TE2-7b_CR"    "REPY1_CR"     "RandI-2"
[36] "RandI-3"      "EnSpm-N2_CR"  "DNA-4-7_CR"   "REP4_CR"      "LTR1_CR"
[41] "TE2-3_CR"     "TE1-1_CR"     "P-1_CR"       "DNA-8-1_CR"   "DNA3-1_CR"
[46] "DNAX1_CR"     "DNA-8-4_CR"   "DNA-3-7_CR"   "DIRS"         "EnSpm-N1a_CR"
[51] "RandI-5"      "TOC2"         "TOC1"         "Gulliver"     "EnSpm-N3_CR"
[56] "TE2-1_CR"     "RandI-4A_CR"  "TE2-2_CR"     "TE2-4_CR"     "DNA-1-9_CR"
[61] "DNA-8-5_CR"   "DNA-8-2_CR"   "REPX4_CR"     "RandI-4B_CR"  "REP1_CR"
[66] "NonLTR-3_CR"  "NonLTR-2_CR"  "LTR3_CR"      "TE2-6_CR"     "TCR1_LTR"
[71] "REM1-LTR"     "EnSpm-N1_CR"  "REPX3_CR"     "DNA-8-6_CR"   "DNA-8-3_CR"
[76] "REPX1_CR"     "Harbinger"

export.gff3(annot,"/home/nem34/annot_newtrans.gff3")
#annot<-import.gff3("/home/nem34/annot_newtrans.gff3")

#Get rid of ones to remove
x<-grep("REP",annot$Name)
annot<-annot[-x] #lost 2000
x<-grep("TE2auto_CR",annot$Name)
annot<-annot[-x] #8870 left

#Make file of SoloLTRs
w<-grep("LTR1",annot$Name)
x<-grep("LTR2",annot$Name)
y<-grep("LTR3",annot$Name)
z<-grep("LTR4",annot$Name)

total<-c(w,x,y,z)
annotSoloLTR <- annot[total]
annotSoloLTR$Type <- "SoloLTR"
export.gff3(annotSoloLTR,"/home/nem34/SoloLTR.gff3")

annot<-annot[-total]

#8660 left

x<-grep("DNA-1",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNA-2",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNA-3",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNA3",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNAX1",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNA-4",annot$Name)
annot$Name[x]<-"DNA.NonAut"
x<-grep("DNA-8",annot$Name)
annot$Name[x]<-"TIR.P"
x<-grep("Gypsy",annot$Name)
annot$Name[x]<-"LTR.Gypsy"
x<-grep("Gulliver",annot$Name)
annot$Name[x]<-"TIR.Gulliver"
x<-grep("Harbinger",annot$Name)
annot$Name[x]<-"TIR.Harbinger"
x<-grep("EnSpm",annot$Name)
annot$Name[x]<-"TIR.EnSpm"
x<-grep("DIRS",annot$Name)
annot$Name[x]<-"DIRS"
x<-grep("hAT",annot$Name)
annot$Name[x]<-"TIR.hAT"
x<-grep("Copia",annot$Name)
annot$Name[x]<-"LTR.Copia"
x<-grep("L1",annot$Name)
annot$Name[x]<-"LINE.L1"
x<-grep("SINEX",annot$Name)
annot$Name[x]<-"SINE"
x<-grep("TE1",annot$Name)
annot$Name[x]<-"Un.TE1"
x<-grep("TE2",annot$Name)
annot$Name[x]<-"Un.TE2"
x<-grep("TOC1",annot$Name)
annot$Name[x]<-"LTR.TOC1"
x<-grep("TOC2",annot$Name)
annot$Name[x]<-"TIR.TOC2"
x<-grep("TCR1",annot$Name)
annot$Name[x]<-"TIR.TCR1"
x<-grep("Novosib",annot$Name)
annot$Name[x]<-"TIR.Novosib"
x<-grep("REM1",annot$Name)
annot$Name[x]<-"LTR.REM1"
x<-grep("RTE",annot$Name)
annot$Name[x]<-"LINE.RTE"
x<-grep("RandI",annot$Name)
annot$Name[x]<-"LINE.Dualen/RandI"
x<-grep("P-1",annot$Name)
annot$Name[x]<-"TIR.P"
x<-grep("P-2",annot$Name)
annot$Name[x]<-"TIR.P"
x<-grep("NonLTR",annot$Name)
annot$Name[x]<-"NonLTR"
x<-grep("nonLTR",annot$Name)
annot$Name[x]<-"NonLTR"
x<-grep("Mariner",annot$Name)
annot$Name[x]<-"TIR.Mariner"

#still 8660

unique(annot$Name)
 [1] "NonLTR"            "LINE.L1"           "LTR.Gypsy"
 [4] "LTR.Copia"         "TIR.hAT"           "Un.TE1"
 [7] "LINE.RTE"          "DNA.NonAut"        "Un.TE2"
[10] "TIR.Novosib"       "LINE.Dualen/RandI" "TIR.P"
[13] "TIR.Mariner"       "LTR.REM1"          "SINE"
[16] "TIR.EnSpm"         "DIRS"              "TIR.TOC2"
[19] "LTR.TOC1"              "TIR.Gulliver"      "TIR.TCR1"
[22] "TIR.Harbinger"

table(annot$Name)

             DIRS        DNA.NonAut LINE.Dualen/RandI           LINE.L1
               92               649               299              2027
         LINE.RTE         LTR.Copia         LTR.Gypsy          LTR.REM1
              459               240               861                57
           NonLTR              SINE         TIR.EnSpm      TIR.Gulliver
              916               125               137                34
    TIR.Harbinger           TIR.hAT       TIR.Mariner       TIR.Novosib
                1              1145               370               372
            TIR.P          TIR.TCR1          TIR.TOC2          LTR.TOC1
              176                 9                13                23
           Un.TE1            Un.TE2
              360               295

			  
annot$type <- annot$Name
export.gff3(annot,"/home/nem34/transposons_newcut.gff3")
#annot<-import.gff3("/home/nem34/transposons_newcut.gff3")

#For selectP
c("LTR.Copia","LTR.Gypsy","LTR.TOC1","LTR.REM1","DIRS","LINE.RTE","LINE.L1","LINE.Dualen/RandI","SINE",
	"TIR.TCR1","TIR.TOC2","TIR.hAT","TIR.Novosib","TIR.Gulliver","TIR.Harbinger","TIR.Mariner","TIR.P","TIR.EnSpm",
	"DNA.NonAut","Un.TE1","Un.TE2","NonLTR")                                    
                            

#Newest edit							
annot<-import.gff3("/home/nem34/transposons_newcut.gff3")


#Set Classes
x<-grep("DNA",annot$type)
annot$class[x]<-"DNA"
x<-grep("TIR",annot$type)
annot$class[x]<-"DNA"
x<-grep("LTR",annot$type)
annot$class[x]<-"RET"
x<-grep("DIRS",annot$type)
annot$class[x]<-"RET"
x<-grep("LINE",annot$type)
annot$class[x]<-"RET"
x<-grep("SINE",annot$type)
annot$class[x]<-"RET"
x<-grep("Un.TE1",annot$type)
annot$class[x]<-"Un"
x<-grep("Un.TE2",annot$type)
annot$class[x]<-"Un"

#Set Orders
x<-grep("DNA",annot$type)
annot$order[x]<-"Un"
x<-grep("NonLTR",annot$type)
annot$order[x]<-"Un"
x<-grep("Un",annot$type)
annot$order[x]<-"Un"
x<-grep("DNA",annot$type)
annot$order[x]<-"Un"
x<-grep("LINE",annot$type)
annot$order[x]<-"LINE"
x<-grep("TIR",annot$type)
annot$order[x]<-"TIR"
x<-grep("DIRS",annot$type)
annot$order[x]<-"DIRS"
x<-grep("SINE",annot$type)
annot$order[x]<-"SINE"
x<-grep("LTR.Copia",annot$type)
annot$order[x]<-"LTR"
x<-grep("LTR.Gypsy",annot$type)
annot$order[x]<-"LTR"
x<-grep("LTR.TOC1",annot$type)
annot$order[x]<-"LTR"
x<-grep("LTR.REM1",annot$type)
annot$order[x]<-"LTR"

#All 8660 annotated

export.gff3(annot,"/home/nem34/transposons_newcut.gff3")
#import.gff3(annot,"/home/nem34/transposons_newcut.gff3")


  





