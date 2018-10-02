#!/usr/bin/python
"""
19th of April 2012
Script to take the gff file containing DE loci
by Sebastian and align several sRNA libraries 
against them.The alignments will then be used
to perform a phasing analysis on each specific 
loci.

author:bacms2
"""
from time import ctime
from sys import argv
from FileOperator import RNAsignature
from Methods import Methods
if __name__=='__main__':
    print ctime()
    #Get the file containing the loci
    zFile = argv[1]
    #Open the file
    handle = open(zFile)
    #Process each line into a dictionary
    lsloci = []
    for zline in handle:
        if zline[0]!='#':
            lsline = zline.strip('\n').split('\t')
            zchr,empty,empty,nstart,nend,empty,empty,empty,zname = lsline
            lsloci.append(('%s,%s,%s' % (zchr,nstart,nend),[zchr[-1],nstart,nend]))
            #dtloci['%s,%s,%s' % (zchr,nstart,nend)] = [zchr,nstart,nend]
    handle.close()
    
    #Define the list of sRNA libraries to be used
    lslibs = ['SL258','SL259','SL302','SL303','SL362','SL363']
    
    #Instantiate the method to calculate the phasing scores
    phasing = Methods()
    
    #Iterate over each library and for each loci store the dtSignatures 
    matrix_results = [[0]*len(lslibs) for i in xrange(len(lsloci))]
    print len(matrix_results)
    loci_pos = 0
    dtsignatures = {}
    lib_pos = 0
    for zlib in lslibs:
        print 'Now processing %s' % zlib
        ls_srnas = []
        #Read the alignment file and store the sRNAs into a list of signatures objects
        handle = open('../datasets/%s.v_Arabidopsis_thaliana_genome-tair9.patman.gff2' % zlib)
        for zline in handle:
            if zline[0]!='#':
                lsline = zline.rstrip('\n').split('\t')
                zchr,empty,empty,nstart,nend,ncount,zstrand,empty,ztag = lsline
                zstrand+='1'
                if zchr!='chloroplast' and zchr!='mitochondria':ls_srnas.append(RNAsignature(zchr,nstart,zstrand,ztag,ncount,len(ztag)))
        #Given the loci and the sRNAs find overlaps between them
        loci_pos = 0#store the position on the list for the loci
        srna_pos = 0#store the position on the list for the sRNAs
        while 1:
            #print loci_pos,srna_pos
            #Ensure both the loci and the sRNAs are on the same chromossome
            while ls_srnas[srna_pos].zChr!=lsloci[loci_pos][1][0]:
                if int(ls_srnas[srna_pos].zChr)<lsloci[loci_pos][1][0]:srna_pos+=1
                else:loci_pos+=1
                
            #First test whether the sRNA end is before the end of the locus
            if (ls_srnas[srna_pos].nCoord+ls_srnas[srna_pos].nSize-1)<=int(lsloci[loci_pos][1][2]):
                #If it is before the end check whether it's start coordinate is inside the locus and save that sRNA
                if ls_srnas[srna_pos].nCoord>=int(lsloci[loci_pos][1][1]):
                    srna = ls_srnas[srna_pos]
                    dtsignatures['%s,%s,%s' %  (srna.zChr,srna.nCoord,srna.nStrand)]=srna
                srna_pos+=1
                     
            else:
                #If the sRNA being consider is after the end of the locus get the next locus
                score = phasing.method1(dtsignatures,24,0,63,105)
                matrix_results[loci_pos][lib_pos]=score[0]
                dtsignatures={}
                loci_pos+=1
                
            #Break the while if no more sRNAs or no more loci
            if loci_pos==len(lsloci) or srna_pos==len(ls_srnas):break
        #Move to the next lib
        lib_pos+=1
    
    handle = open('Sebastian loci with phasing scores.csv','w')
    for zlib in lslibs:handle.write('\t%s'%zlib)
    handle.write('\n')
    for loci_pos in xrange(len(lsloci)):
        handle.write('%s\t' % lsloci[loci_pos][0])
        for lib_pos in xrange(len(lslibs)):
            handle.write('%s\t' % (matrix_results[loci_pos][lib_pos]))
        handle.write('\n')
    handle.close()
    #for key in dt_results['SL258']:print len(dt_results['SL258'][key])
    print ctime()
            
        