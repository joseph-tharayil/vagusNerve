# Apache-2.0

import numpy as np
import pandas as pd

import quantities as pq

from vagusNerve.nerveSetup import *
from vagusNerve.recruitment import *

def PhiWeight(d, current,fascIdx, fascTypes,stimulusDirectory, distribution_params, variance=0,outputfolder=None):

    '''
    For a given current level, for each fiber type calculates the scaling factor for each diameter, based on recruitment, fiber diameter distribution, and number of fibers
    '''


    phiWeight = [ [[],[]], [[],[]] ]

    recruitment = Recruitment(current,d,fascIdx,stimulusDirectory, variance)

    maffProb, meffProb, ueffProb, uaffProb = getFiberTypeFractions(fascIdx, fascTypes,distribution_params)

    if outputfolder is not None:

        np.save(outputfolder + '/fiberTypeFractions_'+str(fascIdx)+'.npy',[maffProb,meffProb])

    numFibersPerFascicle = getFibersPerFascicle(fascIdx,fascTypes,distribution_params)


##### Weight is given by the product of the recruitment curve and the diameter probability curve
    phiWeight[0][0] =  MaffProb(d,maffProb,distribution_params,fascIdx)  * recruitment * numFibersPerFascicle
    phiWeight[0][1] =  MeffProb(d,meffProb,distribution_params,fascIdx)  * recruitment * numFibersPerFascicle

    if outputfolder is not None:
        np.save(outputfolder+'/maffProb'+str(fascIdx)+'.npy',MaffProb(d,maffProb,distribution_params,fascIdx))
        np.save(outputfolder + '/meffProb' + str(fascIdx) + '.npy', MeffProb(d, maffProb, distribution_params, fascIdx))
        np.save(outputfolder+'/numFibersPerFascicle'+str(fascIdx)+'.npy',numFibersPerFascicle)

    return phiWeight, recruitment

def getPhiWeight(d, current,fascIdx,fascTypes, stimulusDirectory, distribution_params, variance=np.array([0]),outputfolder=None):

    '''
    Iterates through current levels and, for each fiber type, returns scaling factor given by product of recruitment, fiber diameter distribution, and number of fibers
    '''

    phiWeight = []
    recruitment = []

    for c in current:

        for v in variance:

            p, rec = PhiWeight(d,c,fascIdx,fascTypes,stimulusDirectory, distribution_params, v,outputfolder)

            phiWeight.append(p)
            recruitment.append(rec)



    phiWeight0 = phiWeight[0][0][0][np.newaxis]
    phiWeight1 = phiWeight[0][0][1][np.newaxis]


    for i in np.arange(1,len(phiWeight)):
        phiWeight0 = np.vstack((phiWeight0,phiWeight[i][0][0][np.newaxis]))
        phiWeight1 = np.vstack((phiWeight1,phiWeight[i][0][1][np.newaxis]))

    if outputfolder is not None:
        np.save(outputfolder+'/recruitment_'+str(fascIdx)+'.npy',recruitment)

    return phiWeight0, phiWeight1
