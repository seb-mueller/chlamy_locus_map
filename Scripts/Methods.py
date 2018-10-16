"""
Python script with implementation of the current methods for phasing detection
author: bacms2
date: 03/08/2011
"""
from phasingRemodel import libraryManager
from FileOperator import RNAsignature
#from mpmath import binomial,hyper
#from numpy import log
import math
#try:
import rpy2.rinterface as rinterface
rinterface.initr()
phyper = rinterface.globalenv.get("phyper")
bolean_rpy2 = True
myparams = {'lower.tail':rinterface.SexpVector([False,],rinterface.LGLSXP),'log.p':rinterface.SexpVector([True,],rinterface.LGLSXP)}
#except ImportError:bolean_rpy2 = False


class Methods():
    def __init__(self):
        pass
    
    """
    def mhyper(self,q,m,n,k,nSize):
        u = min(nSize,n)-1
        A = binomial(m,q)*binomial(n,k-q)
        B = hyper([1,q-k,q-m],[1+q,1-k+n+q], 1)
        C = -binomial(m,1+u)*binomial(n,k-1-u)
        D = hyper([1,1-k+u,1-m+u],[2+u,2-k+n+u], 1)
        return log(float((A*B + C*D) / binomial(m+n, k)))
    """
    
    def phaser(self,dtSigns,nSize,pValue,nMinLength,nMaxGap):
        """
        M.phaser(dtSigns,nSize,pValue,nMinLength,nMaxGap)
        Runs bphase method
        """
        manager = libraryManager(dtSigns,nSize)
        manager.createNodes(nMaxGap)
        #print len(manager.dtSignatures)
        lsValues=manager.computeSignificance(pValue,nMinLength)
        #print lsValues,nResults
        return lsValues
    
    def chen(self,dtSigns,nSize,nCycles):
        """
        M.method2(dtSigns,nSize,nCycles)
        Runs Ho-Ming method return the minimum p-value
        """
        #The window size for searching
        nWindow = nCycles * nSize
        nResults=0
        lsValues = ['0.0','None','None','None']
        for signature in dtSigns.values():
            n=-1
            k=-1
            if signature.nStrand == 1:
                #print 'Signature coordinate:%d'%signature.nCoord
                nAbsCoord = signature.nCoord+nWindow-1
                j=signature.nCoord
                l=signature.nCoord+nSize-3
                for i in xrange(signature.nCoord,nAbsCoord+1):
                    #Calculate n on the sense strand
                    zName = signature.zChr+','+str(i)+','+'1'
                    try:
                        dtSigns[zName]
                        n+=1
                    except KeyError:pass
                    #Calculate k on the same strand
                    if i==j:
                        #print 'i,j',i,j
                        zName = signature.zChr+','+str(i)+','+'1'
                        try:
                            dtSigns[zName]
                            k+=1
                        except KeyError:pass
                        j+=nSize
                    #Calculate n on the anti sense strand
                    if i<=nAbsCoord:
                        zName = signature.zChr+','+str(i-2)+','+'-1'
                        try:
                            dtSigns[zName]
                            n+=1
                        except KeyError:pass
                    #Calculate k on the anti sense strand
                    if i==l and i<=nAbsCoord-2:
                        zName = signature.zChr+','+str(i)+','+'-1'
                        try:
                            dtSigns[zName]
                            k+=1
                        except KeyError:pass
                        l+=nSize    
            elif signature.nStrand==-1:
                nAbsCoord = signature.nCoord-nWindow+1
                j=nAbsCoord+nSize-1
                l=nAbsCoord+2
                for i in xrange(nAbsCoord,signature.nCoord+3):
                    #Calculate n on the sense strand
                    if i<=signature.nCoord:
                        zName = signature.zChr+','+str(i)+','+'-1'
                        try:
                            dtSigns[zName]
                            n+=1
                        except KeyError:pass
                    #Calculate k on the same strand
                    if i==j and i<=signature.nCoord:
                        zName = signature.zChr+','+str(i)+','+'-1'
                        try:
                            dtSigns[zName]
                            k+=1
                        except KeyError:pass
                        j+=nSize
                    #Calculate n on the anti sense strand
                    if i<=signature.nCoord:
                        zName = signature.zChr+','+str(i+2)+','+'1'
                        try:
                            dtSigns[zName]
                            n+=1
                        except KeyError:pass
                    #Calculate k on the anti sense strand
                    if i==l and i<=signature.nCoord+2:
                        zName = signature.zChr+','+str(i)+','+'1'
                        try:
                            dtSigns[zName]
                            k+=1
                        except KeyError:pass
                        l+=nSize
            
            #Calculate p-value from n and k
            phyper_q = rinterface.SexpVector([k-1,], rinterface.INTSXP)
            phyper_m = rinterface.SexpVector([nCycles*2-1,], rinterface.INTSXP)
            phyper_n = rinterface.SexpVector([nWindow*2-(nCycles*2),], rinterface.INTSXP)
            phyper_k = rinterface.SexpVector([n,], rinterface.INTSXP)
            fPvalue =  phyper(phyper_q,phyper_m,phyper_n,phyper_k,**myparams)[0]
        

            if signature.nStrand==1:
                if fPvalue<=float(lsValues[0]):
                    lsValues=[str(fPvalue),str(signature.zChr),str(signature.nCoord),str(signature.nCoord+nWindow)]
            elif signature.nStrand==-1:
                if fPvalue<=float(lsValues[0]):lsValues=[str(fPvalue),str(signature.zChr),str(signature.nCoord-nWindow),str(signature.nCoord)]
            nResults+=1
        return lsValues
    
    
    def brachy(self,dtSigns,nSize,nCycles=10,bFN=False):
        dtSigns1 = dtSigns.copy()
        lsValues = ['0.0','None','None','None']
        amp_factor = 10
        nWindow = nCycles*nSize
        #Move the count from the antisense strand to the sense strand
        for zName,signature in dtSigns1.items():
            if signature.nStrand==-1:
                nCoord = signature.nCoord - (nSize-3)
                #Update counts if position has already matches and create a new entry if not
                try: dtSigns1['%s,%d,1' % (signature.zChr,nCoord)].nCounts+=signature.nCounts
                except KeyError:dtSigns1['%s,%d,1' % (signature.zChr,nCoord)]=RNAsignature(signature.zChr,nCoord,1,nCoord,signature.nCounts,nSize)
                dtSigns1.pop(zName)
        #For every sRNA on the library
        for signature in dtSigns1.itervalues():
            #if bFN and (signature.nCoord-100)%nSize!=0:continue
            window_end = signature.nCoord+nCycles*nSize
            hash_3 = {}#Define an hash to store number of matches in phase
            phased = []#List to store the number of counts in phase
            non_phased = []#List to store the number of counts out of the phase register
            for i in xrange(signature.nCoord,window_end):#Iterate over the positions inside the window
                zName = signature.zChr+','+str(i)+','+'1'#Get the name of the sequence for the position
                try:temp_signature = dtSigns1[zName]#Check if the position has sRNAs
                except KeyError:continue
                #If position divided by the size class has reminder zero is in phase
                if ((temp_signature.nCoord-signature.nCoord)%nSize==0):
                    hash_3[i]=1
                    phased.append(temp_signature.nCounts)
                elif ((temp_signature.nCoord+1-signature.nCoord)%nSize==0):
                    hash_3[i+1]=1
                    phased.append(temp_signature.nCounts)
                elif ((temp_signature.nCoord-1-signature.nCoord)%nSize==0):
                    hash_3[i-1]=1
                    phased.append(temp_signature.nCounts)
                #If position not in phase
                else:non_phased.append(temp_signature.nCounts)
            #print phased,non_phased
            sum_ki=sum(phased)#Get number of reads in phasing
            kp_max = max(phased)#Get position with maximum number of reads
            sum_non_phased = sum(non_phased)#Get number of reads out of phasing phasing
            k_p = len(hash_3)#Get number of phase position occupied
            #print kp_max,sum_ki
            if k_p>=3 and (kp_max/sum_ki)<=0.9:#If window has at least 3 positions occupied calculate phase score
                fPvalue = math.log((1+amp_factor*(sum_ki)/(1+float(sum_non_phased)))**(k_p-2))
                if signature.nStrand==1:
                    if fPvalue>=float(lsValues[0]):lsValues=[str(fPvalue),str(signature.zChr),str(signature.nCoord),str(signature.nCoord+nWindow)]
                elif signature.nStrand==-1:
                    if fPvalue>=float(lsValues[0]):lsValues=[str(fPvalue),str(signature.zChr),str(signature.nCoord-nWindow),str(signature.nCoord)]
            #if pValue>=25.0 and kp_max/sum_ki<=0.9:nResults+=1 # If pValue>Threshold increase number of positive results
        return lsValues
