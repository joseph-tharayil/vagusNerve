# Apache-2.0

import numpy as np
import pandas as pd

import quantities as pq

from vagusNerve.nerveSetup import *
from vagusNerve.recruitment import *

def PhiWeight(d, current,fascIdx, fascTypes,stimulusDirectory, distribution_params, variance=0):

    '''
    For a given current level, for each fiber type calculates the scaling factor for each diameter, based on recruitment, fiber diameter distribution, and number of fibers
    '''


    phiWeight = [ [[],[]], [[],[]] ]

    recruitment = Recruitment(current,d,fascIdx,stimulusDirectory, variance)

    maffProb, meffProb, ueffProb, uaffProb = getFiberTypeFractions(fascIdx, fascTypes, distribution_params)

    numFibersPerFascicle = getFibersPerFascicle(fascIdx,fascTypes)


##### Weight is given by the product of the recruitment curve and the diameter probability curve
    phiWeight[0][0] =  MaffProb(d,maffProb)  * recruitment[0] * numFibersPerFascicle
    phiWeight[0][1] =  MeffProb(d,meffProb)  * recruitment[0] * numFibersPerFascicle


    phiWeight[1][0] =  UaffProb(d,uaffProb)  * recruitment[-1] * numFibersPerFascicle
    phiWeight[1][1] =  UeffProb(d,ueffProb)  * recruitment[-1] * numFibersPerFascicle

    return phiWeight, recruitment

def getPhiWeight(d, current,fascIdx,fascTypes, stimulusDirectory, distribution_params, variance=np.array([0])):

    '''
    Iterates through current levels and, for each fiber type, returns scaling factor given by product of recruitment, fiber diameter distribution, and number of fibers
    '''

    phiWeight = []
    recruitment = []

    for c in current:

        for v in variance:

            p, rec = PhiWeight(d,c,fascIdx,fascTypes,stimulusDirectory, distribution_params, v)

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
