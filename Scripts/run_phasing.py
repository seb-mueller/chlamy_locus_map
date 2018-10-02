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

import csv


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
    lslibs = []
    with open('/projects/nick_mattews/chlamy_locus_map_github/Summary_of_Data.csv') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            if row['InCurrentLociRun'] == "Yes":
                lslibs.append(row['File'])

    #Define the size of sRNAs
    nsize = 21
    
    #Define location of the file containing the loci definition
    zname = '/projects/nick_matthews/phasing/loci_for_phasing.csv'
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
        
    handle = open('/projects/nick_matthews/phasing_results_by_locus_%dnt.tsv' % nsize,'w')
    
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
            
