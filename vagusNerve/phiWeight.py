import numpy as np
import pandas as pd

import quantities as pq

from .utils import *
from .nerveSetup import getFiberTypeFractions
from .phiShape import *
from .recruitment import *

def PhiWeight(d, current,fascIdx, fascTypes,stimulusDirectory):
    
    
    phiWeight = [ [[],[]], [[],[]] ]
    
    recruitment = Recruitment(current,d,fascIdx,stimulusDirectory)
    
    scaling = []
    
    scalingFactors = [1,2]
    
    maffProb, meffProb, ueffProb, uaffProb = getFiberTypeFractions(fascIdx, fascTypes)
    
    numFibersPerFascicle = getFibersPerFascicle(fascIdx,fascTypes)
    
    
##### Weight is given by the product of the recruitment curve and the diameter probability curve
    phiWeight[0][0] =  MaffProb(d,maffProb)  * recruitment[0] * numFibersPerFascicle
    phiWeight[0][1] =  MeffProb(d,meffProb)  * recruitment[0] * numFibersPerFascicle
    
    
    phiWeight[1][0] =  UaffProb(d,uaffProb)  * recruitment[-1] * numFibersPerFascicle
    phiWeight[1][1] =  UeffProb(d,ueffProb)  * recruitment[-1] * numFibersPerFascicle

    return phiWeight,recruitment

def getPhiWeight(d, current,fascIdx,fascTypes, stimulusDirectory):
    
    phiWeight = []
    recruitment = []
    
    for c in current:
        
        p, rec = PhiWeight(d,c,fascIdx,fascTypes,stimulusDirectory)

        phiWeight.append(p)
        recruitment.append(rec)
        

    phiWeight0 = phiWeight[0][0][0][np.newaxis]
    phiWeight1 = phiWeight[0][0][1][np.newaxis]
    phiWeight2 = phiWeight[0][1][0][np.newaxis]
    phiWeight3 = phiWeight[0][1][1][np.newaxis]

    for i in np.arange(1,len(phiWeight)):
        phiWeight0 = np.vstack((phiWeight0,phiWeight[i][0][0][np.newaxis]))
        phiWeight1 = np.vstack((phiWeight1,phiWeight[i][0][1][np.newaxis]))
        phiWeight2 = np.vstack((phiWeight2,phiWeight[i][1][0][np.newaxis]))
        phiWeight3 = np.vstack((phiWeight3,phiWeight[i][1][1][np.newaxis]))
    
    return phiWeight0, phiWeight1, phiWeight2, phiWeight3