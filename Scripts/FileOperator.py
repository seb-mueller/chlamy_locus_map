#!/usr/bin/python
"""
Created on Jan 19, 2010
Implementation of parser for default, gff and patman alignment files
and the structure to store each RNA signature data
@author: bacms2
"""
import os
from sys import exit
class RNAsignature():
    def __init__(self,zChr,nCoord,nStrand,zTag,nCounts,nSize,nEdits=0):
        """
        RNAsignature(self,zChr,nCoord,nStrand,zTag,nCounts,nSize)
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
        self.nSize = nSize
        self.nEdits = nEdits
    
    def __str__(self):
        """
        self.__str__(self)
        Define a human readable string for the object
        """
        return ('chr:%s coord:%d str:%d tag:%s count:%d')%(self.zChr,self.nCoord,self.nStrand,self.zTag,self.nCounts)
    
class alignmentFile():
    def __init__(self,zfInput,nSize=21,zfFormat='default',bIgnoreSize=False):
        """
        Read the alignment files and return a dictionary with all 
        the reads. Each value of the map is an RNAsignature object.
        zfInput ->Path to the input file
        nSize -> Size os sRNAs to search for
        fFormat -> The format of the alignment file
        """
        self.dtSignatures = {}
        self.readFile(zfInput,nSize,zfFormat,bIgnoreSize)
        
        
    def readFile(self,zfInput,nSize,zfFormat,bIgnoreSize = False):
        """
        A.readFile(zfInput,fFormat,nSize)
        The read file function itself, it is called in the instantiation
        of the class and should not be called directly
        """
        zfFormat = zfFormat.lower()
        handle = open(zfInput)
        #default format
        if zfFormat=='default':
            for zLine in handle:
                zLine = zLine.replace(os.linesep,'')
                if zLine == '':continue
                lsPar = zLine.split('\t')
                if bIgnoreSize:#Ignore the size of the sRNA
                    try:self.dtSignatures[lsPar[0]+','+lsPar[1]+','+lsPar[2]].nCounts+=int(lsPar[4])
                    except KeyError:self.dtSignatures[lsPar[0]+','+lsPar[1]+','+lsPar[2]]=RNAsignature(lsPar[0],lsPar[1],lsPar[2],lsPar[3],lsPar[4],nSize)
                elif len(lsPar[3])==nSize:
                    try:self.dtSignatures[lsPar[0]+','+lsPar[1]+','+lsPar[2]].nCounts+=int(lsPar[4])
                    except KeyError:self.dtSignatures[lsPar[0]+','+lsPar[1]+','+lsPar[2]]=RNAsignature(lsPar[0],lsPar[1],lsPar[2],lsPar[3],lsPar[4],nSize)
                #print lsPar[0],lsPar[1],lsPar[2],lsPar[3],lsPar[4],nSize
        #Patman output default format
        elif zfFormat=='patman': 
            #Creates map to convert - and + to -1 and 1
            dtConvert = {'+':'1','-':'-1'}
            nLines = 0 
            for zLine in handle:
                nLines+=1
                zLine = zLine.rstrip(os.linesep)
                #if zLine=='':continue
                lsPar = zLine.split('\t')
                zChr = lsPar[0].split(' ')[0]
                #print 'zChr',zChr
                #Test the start nPos of the read
                
                if lsPar[4]=='+':nPos = lsPar[2]
                elif lsPar[4]=='-':nPos = lsPar[3]
                try:nCounts = int(lsPar[1].split()[-1].split(':')[-1])
                except ValueError:nCounts = 1
                sequence = lsPar[1].split()[0]
                zName = zChr+','+nPos+','+dtConvert[lsPar[4]]
                #print lsPar[0],nPos,dtConvert[lsPar[5]],sequence,int(nCounts),nSize
                if bIgnoreSize:
                    #Ignore the size of the sRNA
                    #If read is redundant check if it is already in the Hash map and update its number
                    try:self.dtSignatures.has_key([zName]).nCounts+=int(nCounts)
                    except KeyError:self.dtSignatures[zName]=RNAsignature(zChr,nPos,dtConvert[lsPar[4]],sequence,int(nCounts),nSize,lsPar[5])
                elif (int(lsPar[3])-int(lsPar[2])+1)==nSize:#Check if the read has the correct size
                    try:self.dtSignatures[zName].nCounts+=int(nCounts)
                    except KeyError:self.dtSignatures[zName]=RNAsignature(zChr,nPos,dtConvert[lsPar[4]],sequence,int(nCounts),nSize,lsPar[5])
        elif zfFormat == 'out':
            #Creates map to convert - and + to -1 and 1
            dtConvert = {'+':'1','-':'-1'}
            for zLine in handle:
                zLine = zLine.replace(os.linesep,'')
                if zLine=='':continue
                lsPar = zLine.split('\t')
                #If the patman file has the count number as a separate column
                if lsPar[5]=='+':nPos = lsPar[3]
                elif lsPar[5]=='-':nPos = lsPar[4]
                nCounts = int(lsPar[2])
                sequence = lsPar[1]
                
                if bIgnoreSize:#Ignore the size of the sRNA
                    try:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[5]]].nCounts+=int(nCounts)
                    except KeyError:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[5]]]=RNAsignature(lsPar[0],nPos,dtConvert[lsPar[5]],sequence,int(nCounts),nSize)
                elif (int(lsPar[4])-int(lsPar[3])+1)==nSize:#Check if the read has the correct size        
                    try:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[5]]].nCounts+=int(nCounts)
                    except KeyError:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[5]]]=RNAsignature(lsPar[0],nPos,dtConvert[lsPar[5]],sequence,int(nCounts),nSize)
                    
        #patman gff2 output format
        elif zfFormat=='gff2':
            #Creates map to convert - and + to -1 and 1
            dtConvert = {'+':'1','-':'-1'}
            for zLine in handle:
                zLine = zLine.replace(os.linesep,'')
                if zLine=='':continue
                lsPar = zLine.split('\t')
                if lsPar[6]=='+':nPos = lsPar[3]
                elif lsPar[6]=='-':nPos = lsPar[4]
                if bIgnoreSize:
                    try:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[6]]].nCounts+=float(lsPar[5])
                    except KeyError:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[6]]]=RNAsignature(lsPar[0],nPos,dtConvert[lsPar[6]],lsPar[8].split()[-1],lsPar[5],nSize)
                elif int(lsPar[4])-int(lsPar[3])+1==nSize:#Check if the read has the correct size
                    try:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[6]]].nCounts+=float(lsPar[5])
                    except KeyError:self.dtSignatures[lsPar[0]+','+nPos+','+dtConvert[lsPar[6]]]=RNAsignature(lsPar[0],nPos,dtConvert[lsPar[6]],lsPar[8].split()[-1],lsPar[5],nSize)
        else:
            exit("Format not recognized\nPlease refer to the documentations for more details")
        handle.close()#Close File
        return