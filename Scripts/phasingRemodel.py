#!/usr/bin/env python
"""
alpha v0.1 23 Apr 2010
Python Script to perform search of Phased sRNAs
It includes the Node,Arc,Graph, manager classes
and the main code.
@author: bacms2
Changes include:
- The concept of graph to build the loci has been abandoned
- Every position in the genome is considered a node
- No more seed nodes are used
"""
from FileOperator import alignmentFile
import numpy
import time
import sys
import os
from matplotlib import use
use('Agg')
from matplotlib import pyplot
pyplot.ioff() 
import pylab
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Patch
from collections import defaultdict
from math import log
try:from mpmath import binomial,hyper
except ImportError:pass
try: import psyco
except ImportError:pass
try: psyco.full()
except NameError: pass
try:
    import rpy2.rinterface as rinterface
    bolean_rpy2 = True
except ImportError:bolean_rpy2 = False

if bolean_rpy2: #Initialize the R interface if present
    rinterface.initr()
    phyper = rinterface.globalenv.get("phyper")
    SexpVector = rinterface.SexpVector
    myparams = {'lower.tail':SexpVector([False,],rinterface.LGLSXP),'log.p': SexpVector([True,],rinterface.LGLSXP)}
    dtPhyper = {}

class Node:
    def __init__(self,zName,nCounts,zTags,nCoord):
        """
        self.__init__(self,zName,nCounts,zTag,nCoord)
        Creates a new instance of Node.
        A node correspond to position in the genome
        zName <- Chromossome name + Coordinate 
        nCounts <- Tuple of integers corresponding to the count numbers of both strands
        zTags <- Tuple of strings corresponding to the tag sequences of both strands
        nCoord <- Coordinate of the node
        """
        self.zName = zName 
        self.nCounts = nCounts
        self.lsTags = zTags
        self.nCoord = int(nCoord)
    
    def getPar(self):
        """
        self.getPar(self)
        Return the name and coordinate of the signature used to create the node
        """
        lsAux = self.zName.split(',')
        return (lsAux[0],self.nCoord)
    
    def __str__(self):
        """
        self.__str__(self)
        Changes the default representation of the object to a human readable form.
        """
        return ('Name:%s Coord:%d Counts:%d|%d Tags:%s|%s\n' % (self.zName,self.nCoord,self.nCounts[0],self.nCounts[1],self.lsTags[0],self.lsTags[1]))
        
class Graph():
    def __init__(self,dtSignatures):
        """
        self.__init__(self,dtSignatures)
        Creates a new instance of Graph
        dtSignatures <- Dictionary containing all the sRNA signatures in the library
        """
        self.dtNodes = {} #store all the nodes in the graph
        self.dtSignatures = dtSignatures
        self.nSize = 0 #Size of the graph
    
    def addNode(self,zName):
        """
        self.addNode(self,zName)
        add a Node to the graph
        zName <- The name of the node, its normally obtained by joining the Chromossome
        name with its the coordinate on the genome
        eg. chr1_11210 
        """
        try:
            #If node already exists return
            self.dtNodes[zName]
            return
        except KeyError:pass
        #Get the parameters of the node
        nCoord = zName.split(',')[1]
        #If the coordinates are outside the chromosome return
        if nCoord<=0:return
        #else create the node and append it to self.dtNodes
        #First get the sense signature
        nCounts1,zTag1 = (0.0,0.0)
        try:
            nCounts1 = self.dtSignatures[zName+',1'].nCounts
            zTag1 = self.dtSignatures[zName+',1'].zTag
        except KeyError:pass
        #Repeat the same for the opposite strand
        nCounts2,zTag2 = (0.0,0.0)
        try:
            nCounts2 = self.dtSignatures[zName+',-1'].nCounts
            zTag2 = self.dtSignatures[zName+',-1'].zTag
        except KeyError:pass
        self.dtNodes[zName]=Node(zName,(nCounts1,nCounts2),(zTag1,zTag2),nCoord)
        #Increase the number of nodes in the graph
        self.nSize+=1

    def card(self):
        """
        self.card(self)
        return the number of nodes in the graph
        """
        return self.nSize
    
    def displayGraph(self):
        """
        self.displayGraph(self)
        displays a string representation of the graph
        """
        print '\nGraph:'
        for node in self.dtNodes.values():
            print '%s: ' % node.zName, 
            for arc in node.lsOutgoing:
                print '%s(%d), ' % (arc.zTarget,arc.nWeight),
            print
    
class locus:
    def __init__(self,nStartCoord,nEndCoord,nLength,nSize,zChr = 'None'):
        """
        locus(self,nStartCoord,nEndCoord,nLength,nSize)
        Create a new instance of the locus class
        nStartCoord <- The start coordinate of the locus on the genome
        nEndCoord <- The end coordinate of the locus on the genome
        nLength <- The length of the locus
        nSize <- The size class of sRNAs being tested
        zChr <- Optional argument to define the chromossome name
        Note that the length of the locus will always be higher than the distance between
        their start and end coordinates because two rows of zeros are added to the rigth
        and left extremes of the locus for computation purposes
        """
        self.nStartCoord = nStartCoord
        self.nEndCoord = nEndCoord
        self.zChr = zChr
        self.lsSTags = nLength*[0]
        self.aSCounts = numpy.zeros(nLength,dtype='Float64')
        self.lsASTags = nLength*[0]
        self.aASCounts = numpy.zeros(nLength,dtype='Float64')
        self.nLength = nLength
        self.nSize = nSize
        self.pValue = 0.0
        self.tMinIndex = (0,0)
        self.phased_counts = 0
        
    def maskRepeats(self):
        """
        Function to mask the tags that appear more than once on segment.
        """
        #get repeated position on sense strand
        positions_init = (self.aSCounts>0).sum()
        items = defaultdict(list)
        for i,item in enumerate(self.lsSTags):items[item].append(i)
        Indexes = []
        for item, locus in items.iteritems():
            if len(locus) > 1 and item!=0:Indexes.extend(locus)
        self.aSCounts[Indexes]=0
        #if (self.aSCounts>0).sum()<positions_init:print self.zChr, self.nStartCoord,(self.aSCounts>0).sum(),positions_init
        #get repeated position on anti sense strand  
        for key in items.iterkeys():items[key] = []
        positions_init = (self.aASCounts>0).sum()
        for i,item in enumerate(self.lsASTags):items[item].append(i)
        Indexes =[]
        for item, locus in items.iteritems():
            if len(locus) > 1 and item!=0:Indexes.extend(locus)
        #self.aASCounts[Indexes[0]]=0
        self.aASCounts[Indexes]=0
        #if (self.aASCounts>0).sum()<positions_init:print self.zChr, self.nStartCoord,(self.aASCounts>0).sum(),positions_init
        return
    
    def chen(self,lindex,rindex,nMinLength,Indexes=True):
        """
        self.hypergeometric(self,lindex,rindex)
        Performs the hypergeometric test for the locus between [lindex,rindex[
        lindex<-Integer indicating the left position on the array
        rindex <-Integer indicating the right position on the array
        """
        if Indexes:
            #Create the subarray to test
            aASense = self.aASCounts[lindex-2:rindex-2]
            #Test to make sure no mistake is happening with the indexes
            if lindex<0 or rindex<0: raise Exception('Indexes cannot be lower than zero')
            if len(aASense)!= len(self.aSCounts[lindex:rindex]):raise Exception('Sense and Anti sense have different lengths')
            #Stack both strands
            aLocus = numpy.hstack([self.aSCounts[lindex:rindex],aASense[::-1]])
        else: aLocus = numpy.hstack([lindex,rindex[::-1]])
        length = len(aLocus)
        #Get the indexes of Phased positions
        index_nSize=numpy.arange(0,length,self.nSize)
        #Caculate the number of phase positions
        m = length/self.nSize
        #Calculate the number of out of phase positions
        n = (length-1)-(length/self.nSize-1)
        #Calculate the number of phased positions occupied
        q =(aLocus[index_nSize]>0).sum()
        #Calculate the number of positions occupied in the locus
        k = (aLocus>0).sum()
        fPvalue = self.r_phyper(q-1, m, n, k)
        return fPvalue
        
    def bphase(self,lindex,rindex,nMinLength,Indexes=True,loci_extension = False,nCycles=10):
        """
        self.hypergeometric(self,lindex,rindex)
        Performs the hypergeometric test for the locus between [lindex,rindex[ for
        all the count number in phase positions and return the minimum p-Value
        lindex<-Integer indicating the left position on the array
        rindex <-Integer indicating the right position on the array
        """
        if Indexes:
            #Create the subarray to test
            aASense = self.aASCounts[lindex-2:rindex-2]
            #Test to make sure no mistake is happening with the indexes
            if lindex<0 or rindex<0: raise Exception('Indexes cannot be lower than zero')
            if len(aASense)!= len(self.aSCounts[lindex:rindex]):raise Exception('Sense and Anti sense have different lengths')
            #Stack both strands
            if loci_extension:
                extension_length = nCycles*self.nSize - len(self.aSCounts[lindex:rindex])
                if extension_length >0:aLocus = numpy.hstack([self.aSCounts[lindex:rindex],numpy.zeros(2*extension_length),aASense[::-1]])
                else:aLocus = numpy.hstack([self.aSCounts[lindex:rindex],aASense[::-1]])
            else:aLocus = numpy.hstack([self.aSCounts[lindex:rindex],aASense[::-1]])
            
        else: aLocus = numpy.hstack([lindex,rindex[::-1]])
        length = len(aLocus)
        #Get the indexes of Phased positions
        index_nSize=numpy.arange(0,length,self.nSize)
        #Get and array with only the Phased Positions
        aAux =aLocus[index_nSize]
        #Get the unique values in the array of Phased counts
        #lsPhasedValues = numpy.unique1d(aAux[aAux>0]) #Old versions of numpy
        lsPhasedValues = numpy.unique(aAux[aAux>0])

        #Caculate the number of phase and non phase positions in the locus
        m = length/self.nSize
        n = (length-1)-(length/self.nSize-1) 
        #Do a standard call of lphyper to increase speed
        if bolean_rpy2:r_phyper = self.r_phyper
        #Remove the zero values from aLocus to reduce its size
        aLocus = aLocus[aLocus>0]
        nCounts = aLocus.sum()
        #if (aAux>0).sum()<(nMinLength/self.nSize):return 0.0
        #Initialize the p-Value to zero
        fPvalue = 0
        q,k=0,0
        for r in lsPhasedValues:
            #Calculate matches in Phased positions
            q = (aAux>=r).sum()
            #Calculate matches in non Phased positions and reduce the array to values higher than r
            aLocus = aLocus[aLocus>=r]
            k = len(aLocus)
            #Get the p-value for the hypergeometric distribution
            if bolean_rpy2:value = self.r_phyper(q-1, m, n, k)
            else: value = self.mhyper2(q,m,n,k) 
            fPvalue = min(fPvalue,value)
        #if fPvalue<=-10.0:print 'q',q,'m',m,'n',n,'k',k,'p-value',fPvalue
        #Return p-Value
        return fPvalue,nCounts
    
    def brachy(self,lindex,rindex,nMinLength,Indexes=True):
        """
        self.hypergeometric(self,lindex,rindex)
        Performs the hypergeometric test for the locus between [lindex,rindex[
        lindex<-Integer indicating the left position on the array
        rindex <-Integer indicating the right position on the array
        """
        if Indexes:
            #Create the subarray to test
            aASense = self.aASCounts[lindex-2:rindex-2]
            #Test to make sure no mistake is happening with the indexes
            if lindex<0 or rindex<0: raise Exception('Indexes cannot be lower than zero')
            if len(aASense)!= len(self.aSCounts[lindex:rindex]):raise Exception('Sense and Anti sense have different lengths')
            #Stack both strands
            aLocus = self.aSCounts[lindex:rindex]+aASense[::-1]
        
        #Get the length of the locus
        length = len(aLocus)
        #Get indexes of the phased positions
        index_phased=numpy.arange(0,length,self.nSize) #normal phase
        index_phased_plus1 =numpy.arange(1,length,self.nSize) #phased positions +1
        index_phased_minus1 = numpy.arange(self.nSize-1,length,self.nSize) #phased positions -1 
        #index_phased.sort()
        
        #Get number of phase position occupied
        phased = aLocus[index_phased]>0
        phased_plus1 = aLocus[index_phased_plus1]>0
        phased_minus1 = aLocus[index_phased_minus1]>0
        n = (phased+phased_plus1+phased_minus1).sum()
        
        #Get number of sRNAs in phase
        p =(aLocus[index_phased[aLocus[index_phased]>0]]).sum()
        p+=(aLocus[index_phased_plus1[aLocus[index_phased_plus1]>0]]).sum()
        p+=(aLocus[index_phased_minus1[aLocus[index_phased_minus1]>0]]).sum()
        
        #Get number of reads out of phasing
        u = (aLocus[aLocus>0]).sum()-p
        
        #Calculate phasing score
        if n>3:fPvalue = log(pow(1+10*(p/(1+u)),n-2))
        else: fPvalue=0
        return fPvalue
    
    def mhyper(self,q,m,n,k):
        """
        loc.mhyper(q,m,n,k)
        Calculate p-value using binomial function from mpmath
        """
        return numpy.log(float(sum(binomial(m,i)*binomial(n,k-i)/binomial(n+m,k) for i in xrange(q,self.nSize))))
    
    def mhyper2(self,q,m,n,k):
        u = min(self.nSize,n)-1
        A = binomial(m,q)*binomial(n,k-q)
        B = hyper([1,q-k,q-m],[1+q,1-k+n+q], 1)
        C = -binomial(m,1+u)*binomial(n,k-1-u)
        D = hyper([1,1-k+u,1-m+u],[2+u,2-k+n+u], 1)
        return numpy.log(float((A*B + C*D) / binomial(m+n, k)))
    
    def r_phyper(self,q,m,n,k):
        """
        self.phyper(self,q,m,n,k)
        Calculate p-value using R function phyper from rpy2 low-level
        interface. 
        "R Documentation
        phyper(q, m, n, k, lower.tail = TRUE, log.p = FALSE)
        q: vector of quantiles representing the number of white balls
            drawn without replacement from an urn which contains both
            black and white balls.
        m: the number of white balls in the urn.
        n: the number of black balls in the urn.
        k: the number of balls drawn from the urn.
        log.p: logical; if TRUE, probabilities p are given as log(p).
        lower.tail: logical; if TRUE (default), probabilities are P[X <= x],
            otherwise, P[X > x].
        "
        """
        phyper_q = SexpVector([q,], rinterface.INTSXP)
        phyper_m = SexpVector([m,], rinterface.INTSXP)
        phyper_n = SexpVector([n,], rinterface.INTSXP)
        phyper_k = SexpVector([k,], rinterface.INTSXP)
        #print q,m,n,k,phyper(phyper_q,phyper_m,phyper_n,phyper_k,**myparams)[0]
        return phyper(phyper_q,phyper_m,phyper_n,phyper_k,**myparams)[0]
        
    def plot(self,display=False,whole=True):
        """
        self.plot(self,display=False,whole=True)
        all <- Flag can be used to plot the whole locus or just the most significant part
        Plot the number of signatures in the sense and anti sense 
        strand for the current locus.
        """
        nPos = True
        #Create a figure object
        fig = pylab.figure(figsize=(16,9),dpi=160)
        for aStrand in [self.aSCounts,self.aASCounts]:        
            #Define Patch objects for the legend
            non_phase_pos = Patch(edgecolor='b', facecolor='b')        
            phase_pos = Patch(edgecolor='r', facecolor='r')
            ylabel = 'Number of sRNAs'
            #Plot strand
            if whole:#Plot the whole strand
                #Define the title to be shown on the plot
                pylab.suptitle('%s:%d..%d with log p-value = %.2f' % (self.zChr,self.nStartCoord,self.nEndCoord,self.pValue),fontsize=18,x=0.39)
                #Remove the extra elements from the sense strand
                aStrand = aStrand[self.nSize:-self.nSize]
                #get the first position in phase
                nFirstPos = self.tMinIndex[0]-self.tMinIndex[0]/self.nSize*self.nSize
                if not nPos: nFirstPos=nFirstPos-3
                #Get the positions that should contain phased sRNAs
                phasedIndex = numpy.arange(nFirstPos,len(aStrand),self.nSize)
            else: #Plot only the most significant region
                #Define the title to be shown on the plot
                pylab.suptitle('Most significant region on phased sRNA loci in %s:%d..%d with log p-value = %.2f' % (self.zChr,self.nStartCoord+self.tMinIndex[0]-self.nSize,self.nStartCoord+self.tMinIndex[1]-self.nSize,self.pValue),fontsize=12,x=0.45)
                #Remove the extra elements from the sense strand
                aStrand = aStrand[self.tMinIndex[0]:self.tMinIndex[1]]
                #Get the positions that should contain phased sRNAs
                if nPos:phasedIndex = numpy.arange(0,len(aStrand),self.nSize)
                else: phasedIndex = numpy.arange(self.nSize-3,len(aStrand),self.nSize)
            
            #If the values are large plot in log scale
            try:
                if max(aStrand)>=1000: 
                    ylabel = 'log(Number of sRNAs)'
                    aStrand = numpy.log(aStrand)
            except ValueError: pass #May not contain any significant phased region
            
            #Plot title and legend for the figure
            pylab.figlegend([phase_pos,non_phase_pos],['Phased positions','Out of Phase positions'],'upper right',prop={'size':8})
            #Produce the actual plot
            ##Divide the screen in two
            if nPos:pylab.subplot(211)
            else:pylab.subplot(212)
            #Adjust the margins
            pylab.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.92, hspace=0.16)
            #Display the labels on the axis 
            pylab.xlabel('Relative position',fontsize= 14)
            pylab.ylabel(ylabel,fontsize= 14)
            #Plot the bars corrresponding to whole the locus
            pylab.bar(range(len(aStrand)),aStrand,color='b',width=1.0,edgecolor='b')
            #Plot the phased positions bars
            pylab.bar(phasedIndex,aStrand[phasedIndex],color='r',width=2,edgecolor='r',alpha=1,align='center')
            try:yaxis_pos = max(aStrand)/3.0*2 #Compute the position where to place the markers of phased positions
            except ValueError:yaxis_pos = 0
            try:pylab.scatter(phasedIndex, [yaxis_pos]*len(phasedIndex), marker='d', color='white',edgecolor='black')
            except ValueError:pass
            try:#The sense strand may not have any sRNA in phase
                phasedIndexOccupied = phasedIndex[aStrand[phasedIndex]>0]
                pylab.scatter(phasedIndexOccupied, [yaxis_pos]*len(phasedIndexOccupied),s=40,marker='d',color='r',edgecolor='black')
            except ValueError:pass 
            try:pylab.axis([0,len(aStrand),0,max(aStrand)])
            except ValueError:pass
            nPos=False
        if display:
            pylab.show()
            pylab.close()
            return
        return fig

class libraryManager:
    def __init__(self,dtSignatures,nSize):
        """
        libraryManager(self,dtSignatures,nSize)
        dtSignatures <- Dictionary containing signature objects
        nSize <- Size class of sRNAs
        Create a new instance of the libraryManager
        """
        self.graph = Graph(dtSignatures) #Instantiate the graph
        self.lsResults = [] #Will be used to store the Strong Connected Components
        self.nSize = nSize
        self.dtSignatures = dtSignatures
        self.dtlocus = {} #store the results to be used in iPython 
        self.lsPvalues = [] #store the p-values to be used for debug
        
    def createNodes(self,nMaxGap):
        """
        self.createNodes(self,nMaxGap)
        nMaxGap <- Maximum gap allowed between two signatures
        Creates all the nodes based on the signatures on the dtSignatures.
        A node is defined by its chromosomal location and contains information
        about both strands. For the extension step only the presence or absence
        of the node is taken in consideration.
        """
        #Create a temporary dictionary to store the nodes per chromosome 
        dtNodes = {}
        #Create the nodes for all the signatures
        #print 'Number of different sRNA signatures:%d' % len(self.dtSignatures)
        for signature in self.dtSignatures.values():
            zName = '%s,%d' % (signature.zChr,signature.nCoord)
            #Add the node to the graph
            self.graph.addNode(zName)
            #Add the node to the dictionary of its chromossome
            try: dtNodes[signature.zChr].append(signature.nCoord)
            except KeyError:dtNodes[signature.zChr]=[signature.nCoord]
        #Get all the loci
        for zChr in dtNodes:
            scc = []#Define a list to store a locus
            #Sort nodes by their chromosome and then by their coordinate
            dtNodes[zChr] = dict.fromkeys(dtNodes[zChr]).keys()
            dtNodes[zChr].sort()
            #print 'nodes:',dtNodes[zChr]
            #Create the loci
            for i in xrange(0,len(dtNodes[zChr])-1):
                zName = '%s,%d' % (zChr,dtNodes[zChr][i])
                nGap = dtNodes[zChr][i+1]-int(dtNodes[zChr][i])+1
                scc.append(zName)
                if nGap>nMaxGap:#If the gap between two nodes is bigger than the nMaxGap store the current locus and start a new one
                    if len(scc)>1:self.lsResults.append(scc[:])
                    scc=[]
            #Append the last scc since there is no node after to test the distance
            zName = '%s,%d' % (zChr,dtNodes[zChr][len(dtNodes[zChr])-1])
            scc.append(zName)
            if len(scc)>1:self.lsResults.append(scc[:])
        del dtNodes#Call garbage collection for dtNodes to keep a low memory footprint
        return
    
    def buildLocus(self,scc,nMinLength=0,nCycles = 10):
        """
        self.buildLocus(self,scc)
        scc <- List containing the name of the nodes to create a Locus
        Create and return a Locus object, based on the list of Nodes.
        """
        #Calculate the lenght of the region with sRNAs span on it
        nLength = int(scc[-1].split(',')[1])-int(scc[0].split(',')[1])+1
        if nMinLength!=0 and nLength<nMinLength:nLength += nCycles *  self.nSize - nLength
        #Create the locus object, it's 
        loc = locus(int(scc[0].split(',')[1]),int(scc[-1].split(',')[1]),nLength+2*self.nSize,self.nSize)
        loc.zChr,nCoord = scc[0].split(',')
        nCoord = int(nCoord)
        for zName in scc:
            node = self.graph.dtNodes[zName]
            loc.aSCounts[node.nCoord-nCoord+self.nSize],loc.aASCounts[node.nCoord-nCoord+self.nSize] = node.nCounts
            loc.lsSTags[node.nCoord-nCoord+self.nSize],loc.lsASTags[node.nCoord-nCoord+self.nSize] = node.lsTags
        loc.maskRepeats()
        return loc
        
    def computeSignificance_fixed_lenght(self,fThreshold,nMinLength,zfOutput=None,bPdf=False,method='bphase',extension=False):
        #If an output location is specified create the files for them
        if zfOutput:handle = open(zfOutput+'.txt','w')
        if bPdf:pdf = PdfPages(zfOutput+'.pdf')
        lsValues = ['0.0','None','None','None','0']
        if len(self.lsResults)==0:return lsValues
        for scc in self.lsResults:
            #If the segment has less than 2 nodes break
            if len(scc)<=2:continue 
            #Build the locus object by calling the buildLocus method
            loc = self.buildLocus(scc,nMinLength)
            self.dtlocus[scc[0]]=loc
            aMin = 0.0
            tMinIndex = (-1,-1)
            phased_counts = 0
            bPrint = False #Creates a boolean variable to flag if the locus is significant or not
            for i in xrange(2,loc.nLength-self.nSize+2):
                ##If the first node has no reads on phased positions ignore this walk
                #if loc.aSCounts[i]==0 and i+self.nSize-3<loc.nLength and loc.aASCounts[i+self.nSize-3]==0:continue
                
                for j in xrange(nMinLength,loc.nLength,self.nSize):
                    if j-i+1<nMinLength:continue
                    
                    if loc.aSCounts[j-self.nSize]==0 and loc.aASCounts[j-3]==0:continue 
                    if method=='bphase':value,nCounts = loc.bphase(i,j,nMinLength,loci_extension=extension)
                    elif method == 'chen' : value = loc.chen(i,j,nMinLength)
                    elif method == 'brachy':value = loc.brachy(i,j,nMinLength)
                    if method=='bphase' or method == 'chen':
                        if value<aMin:
                            aMin=value
                            tMinIndex = (i,j)
                            phased_counts = nCounts
                    elif value>aMin:
                        aMin=value
                        tMinIndex = (i,j)
            loc.tMinIndex = tMinIndex
            loc.pValue = aMin
            loc.phased_counts = phased_counts
            #If the p-value is higher than the threshold save it on the list
            if method=='bphase' or method == 'chen':
                if loc.pValue<=float(lsValues[0]):
                    nCounts = loc.phased_counts
                    zChr,nCoord= scc[0].split(',')
                    lsValues = [str(loc.pValue),zChr,str(loc.tMinIndex[0]+int(nCoord)-int(self.nSize)),str(loc.tMinIndex[1]+int(nCoord)-int(self.nSize)),nCounts]
                    #lsValues=[str(loc.pValue),str(loc.zChr),str(loc.nStartCoord),str(loc.nEndCoord)]
            elif loc.pValue>=float(lsValues[0]):lsValues=[str(loc.pValue),str(loc.zChr),str(loc.nStartCoord),str(loc.nEndCoord)]
            
            if method=='bphase' or method == 'chen':
                if aMin<=fThreshold:bPrint=True
            elif  aMin>=fThreshold:bPrint=True
            if bPrint and zfOutput:
                bPrint = False
                self.writeOutput(handle,scc,loc)
                if bPdf:
                    #Plot the all locus
                    fig = loc.plot()
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
                    #Plot the phased bit
                    fig = loc.plot(whole=False)
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
        if zfOutput:handle.close()
        if bPdf:pdf.close()
        return lsValues
        
        
    def computeSignificance(self,fThreshold,nMinLength,zfOutput=None,bPdf=False,method='bphase'):
        """
        M.computeSignificance(self,fThreshold,nMinLength,zfOutput,bWrite=False)
        Function to compute hypergeometric tests for all the possible paths 
        connecting the nodes.
        """
        #If an output location is specified create the files for them
        if zfOutput:handle = open(zfOutput+'.txt','w')
        if bPdf:pdf = PdfPages(zfOutput+'.pdf')
        lsValues = ['0.0','None','None','None','0']
        
        if len(self.lsResults)==0:return lsValues
        for scc in self.lsResults:
            #If the segment has less than 2 nodes break
            if len(scc)<=2:continue 
            
            #Build the locus object by calling the buildLocus method
            loc = self.buildLocus(scc)
            
            #Control test to ensure that the scc does not has more sRNAs than 
            if len(scc)>loc.nLength:sys.exit('Error length of scc cannot be higher than loc.nLength')
            
            #If locus is smaller than the minimum length ignore it, mainly for speed purposes 
            if loc.nLength<nMinLength:continue
            
            #Store the locus for posterior analysis under an ipython session
            self.dtlocus[scc[0]]=loc
            #Initiate minimum p-value at 0.0 and phased index at -1,-1
            aMin = 0.0
            tMinIndex = (-1,-1)
            
            #Creates a boolean variable to flag if the locus is significant or not
            bPrint = False 
            
            #Start iterating over each position on the locus as a potential start site for phasing
            for i in xrange(2,loc.nLength-self.nSize+2):
                #If the first node has no reads on phased positions ignore this walk
                if loc.aSCounts[i]==0 and i+self.nSize-3<loc.nLength and loc.aASCounts[i+self.nSize-3]==0:continue
                #Start iterate over potential end positions for phasing
                for j in xrange(i+self.nSize,loc.nLength,self.nSize):
                    #if the locus length is smaller than the minimum size for the locus ignore it
                    if j-i+1<nMinLength:continue
                    #If the new locus as no sRNAs in the end positions ignore this walk 
                    if loc.aSCounts[j-self.nSize]==0 and loc.aASCounts[j-3]==0:continue
                    #Ignore segment if one strand has no counts matching to it
                    #if loc.aSCounts[i:j].sum()==0 or loc.aASCounts[i-2:j-2].sum()==0:continue
                    if method=='bphase':value,nCounts = loc.bphase(i,j,nMinLength)
                    elif method == 'chen' : value = loc.chen(i,j,nMinLength)
                    elif method == 'brachy':value = loc.brachy(i,j,nMinLength)
                    if method=='bphase' or method == 'chen':
                        if value<aMin:
                            aMin=value
                            tMinIndex = (i,j)
                            
                    elif value>aMin:
                        aMin=value
                        tMinIndex = (i,j)
            loc.tMinIndex = tMinIndex
            loc.pValue = aMin
            #If the p-value is higher than the threshold save it on the list
            if method=='bphase' or method == 'chen':
                if loc.pValue<=float(lsValues[0]):
                    #lsValues=[str(loc.pValue),str(loc.zChr),str(loc.nStartCoord),str(loc.nEndCoord)]
                    zChr,nCoord= scc[0].split(',')
                    nCounts = loc.phased_counts
                    lsValues = [str(loc.pValue),str(loc.zChr),str(loc.tMinIndex[0]+int(nCoord)-int(self.nSize)),str(loc.tMinIndex[1]+int(nCoord)-int(self.nSize)),nCounts]
            elif loc.pValue>=float(lsValues[0]):lsValues=[str(loc.pValue),str(loc.zChr),str(loc.nStartCoord),str(loc.nEndCoord)]
            
            if method=='bphase' or method == 'chen':
                if aMin<=fThreshold:bPrint=True
            elif  aMin>=fThreshold:bPrint=True
            if bPrint and zfOutput:
                bPrint = False
                self.writeOutput(handle,scc,loc)
                if bPdf:
                    #Plot the all locus
                    fig = loc.plot()
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
                    #Plot the phased bit
                    fig = loc.plot(whole=False)
                    pdf.savefig(fig)
                    pylab.close(fig)
                    del(fig)
        if zfOutput:handle.close()
        if bPdf:pdf.close()
        return lsValues
                
    def writeOutput(self,handle,scc,loc):
        """ 
        writeOutput(self,handle,scc,loc)
        Write a matrix with the p-values for the locus and a file 
        with a summary of the results
        """
        zChr,nCoord= scc[0].split(',')
        handle.write('Identified sRNA locus in %s from position %s to %i\n' % (zChr,loc.nStartCoord,loc.nEndCoord))
        handle.write('\tPhased detected from position %s to %i\n' % (loc.tMinIndex[0]+int(nCoord)-int(self.nSize),loc.tMinIndex[1]+int(nCoord)-int(self.nSize)))
        handle.write('\tLog p-Value = %f\n' % loc.pValue)
        handle.write('\tTo plot this locus type: manager.dtlocus[\'%s\'].plot()\n' % scc[0])
        handle.write('\tThe coordinates in the graph will be:%i to %i\n' % (loc.tMinIndex[0]-self.nSize,loc.tMinIndex[1]-int(self.nSize)))
        handle.write('%s\n' % ('#'*74))
        return
        
        
if __name__=='__main__':
    print time.ctime()
    print 'Running Phased locus identifier'
    start = time.time()
    if len(sys.argv)<7:
        sys.exit('Not enough parameters for the run\nPlease refer to the documentation for more information')
    #Input file
    zfInput = sys.argv[1]
    #Format of the alignment file
    fFormat = sys.argv[2]
    #P-value to be used
    nThreshold = float(sys.argv[3])
    #if nThreshold>0:sys.exit("###Error###\nThreshold cannot be higher than 0.0")
    #The individual size of each siRNA
    nSize = int(sys.argv[4])
    if nSize<0:sys.exit("###Error###\nSize for small RNAs has to be higher than 0")
    #Minimum length of locus
    nMinLength = int(sys.argv[5])
    #Maximum gap 
    nMaxGap = int(sys.argv[6])
    
    #Reading the input files
    if not os.path.isfile(zfInput):sys.exit('File not found %s' %  zfInput)
    dtSignatures = alignmentFile(zfInput,nSize,fFormat).dtSignatures
    print 'Reading the input file took:%f' % (time.time()-start)
    #Getting the output information
    zSample = zfInput.split('/')[-1]
    print zSample
    #zSample = zfInput.split('.')[0].split('/')[-1].split('.')[0]
    #Open file to store the clusters
    zfOutput = '%s_%d_%s_phasing' % (zSample,nSize,nThreshold)
    if os.path.isfile(zfOutput):print '###Warning###:Output file exists and will be overwritten'
    #try:handleOutput = open(zfOutput, 'w')
    #except:sys.exit('Invalid file or file in used, Output cannot be saved.\nExecution aborted')
    #handleOutput.close()
    
    #Print all the information to the user
    print
    print 'Parameters to be used will be:'
    print 'Input File->%s' % zfInput
    print 'Output File -> %s' % zfOutput
    print 'Length of Phased siRNAs:%d' % nSize
    print 'P-value:%f'% nThreshold
    print 'Number of distinct sRNAs signatures to be used:%d' % len(dtSignatures)
    print 'Minimum length for locus: %d' % nMinLength
    print 'Maximum Gap length: %d' % nMaxGap
    print 'Output file->%s' % zfOutput
    print 
    #Creating the nodes
    manager = libraryManager(dtSignatures,nSize)
    astart = time.time()
    manager.createNodes(nMaxGap)
    print 'Nodes Created\nCreating nodes took:%fs' % (time.time()-astart)
    print 
    astart = time.time()
    manager.computeSignificance(nThreshold,nMinLength,zfOutput,bPdf=True)
    #manager.computeSignificance_fixed_lenght(nThreshold,nMinLength,zfOutput,bPdf=True,method='bphase')
    print 'Computing significance took:%fs' % (time.time()-astart)
    print 'All finished'
    print 'This run took:%f' % (time.time()-start)
    
    
