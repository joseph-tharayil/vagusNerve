# Apache-2.0

import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

from scipy.stats import rv_histogram

from math import gamma

from scipy.optimize import curve_fit

import quantities as pq

def loadDiameterDistribution(fiberType):

    #### Loads diameter distributions

    if fiberType == 'maff':

        vals = np.loadtxt(r'D:\vagusOptimization\Data\maffvals.csv',delimiter=',')
    elif fiberType == 'meff':

        vals = np.loadtxt(r'D:\vagusOptimization\Data\meffvalsSmooth.csv',delimiter=',')

    elif fiberType == 'uaff':
        vals = np.loadtxt(r'D:\vagusOptimization\Data\uaffvals.csv',delimiter=',')

    elif fiberType == 'ueff':
        vals = np.loadtxt(r'D:\vagusOptimization\Data\ueffvals.csv',delimiter=',')

    else:
        raise ValueError('Invalid fiber type')
    #####

    #### Gets midpoints of histogram bins
    binWidths = np.diff(vals[:,0]) # Width of each diameter bin, in um

    diameterBins = (vals[:-1,0] + vals[1:,0])/2 * 1e-6 # Converts from um to m
    probabilities = (vals[:-1,1] + vals[1:,1])/2 * 0.01 # Converts from percentage to fraction
    #######

    probabilities = adjust_probabilities(probabilities, binWidths)

    return diameterBins, probabilities

def adjust_probabilities(probabilities, binWidths):

    '''
    Ensures that the total fiber diameter probability distribution sums to 1
    '''

    relativeBinWidths = binWidths/binWidths[0]

    probabilities *= relativeBinWidths

    probabilities /= np.sum(probabilities)

    return probabilities


def getFiberTypeArea(scalingFactors):

    ''' Gets the area occupied by each fiber type. If byFascicle==False, gets area over the nerve. Else, gets area over each fascicle'''

    maffD, maffP = loadDiameterDistribution('maff')
    meffD, meffP = loadDiameterDistribution('meff')
    uaffD, uaffP = loadDiameterDistribution('uaff')
    ueffD, ueffP = loadDiameterDistribution('ueff')

    maffArea = scalingFactors[0] * np.sum(maffD**2*maffP)
    meffArea = scalingFactors[1] * np.sum(meffD**2*meffP)
    uaffArea = scalingFactors[2] * np.sum(uaffD**2*uaffP)
    ueffArea = scalingFactors[3] * np.sum(ueffD**2*ueffP)

    return maffArea, meffArea, uaffArea, ueffArea

def getFiberTypeArea_Overall():

    '''
    Calculates area occupied by each fiber type over the entire nerve
    '''

    ### Overall fiber counts in the nerve
    maffcount = 34000
    meffcount = 14800
    ueffcount = 21800
    uaffcount = 315000
    ####

    maffArea, meffArea, uaffArea, ueffArea = getFiberTypeArea([maffcount, meffcount, uaffcount, ueffcount])

    return maffArea, meffArea, uaffArea, ueffArea

def getFiberTypeArea_byFascicle(fascIdx, fascTypes):

    '''
    Calculates area occupied by each fiber type in a given fascicle
    '''

    maffFrac, meffFrac, ueffFrac, uaffFrac = getFiberTypeFractions(fascIdx,fascTypes)

    maffArea, meffArea, uaffArea, ueffArea = getFiberTypeArea([maffFrac, meffFrac, uaffFrac, ueffFrac])

    return maffArea, meffArea, uaffArea, ueffArea

def getAreaScaleFactor(fascicleSizes):

    '''
    Finds ratio of total fascicle area to total area occupied by fibers. Assumes that fiber type has no impact on this ratio
    '''

    maffArea, meffArea, uaffArea, ueffArea = getFiberTypeArea_Overall()

    totalFiberArea = maffArea + meffArea + uaffArea + ueffArea

    diamScaleFactor = (np.sum(fascicleSizes)/totalFiberArea)

    return diamScaleFactor

def getNumFibers(fascicleSizes,fascIdx,fascTypes):

    '''
    Calculates total number of fibers in a given fascicle
    '''

    diamScaleFactor = getAreaScaleFactor(fascicleSizes)

    maffArea, meffArea, uaffArea, ueffArea = getFiberTypeArea_byFascicle(fascIdx,fascTypes)

    fascicleNumber = fascicleSizes[fascIdx] / (diamScaleFactor * (maffArea + meffArea + uaffArea + ueffArea))

#    fascicleNumber = 10000

    return fascicleNumber


def sampleFractionHistogram(ColorX,ColorY,Colors,side,rng):

    ##### In this function, we generate a histogram of per-fascicle fiber type fractions from the paper, and sample from it to produce the fraction for the simulated fascicle

    Interp = interp1d(ColorX,ColorY,bounds_error=False,fill_value = 'extrapolate') # Interpolation object for the color scale bar
    FracList = Interp(Colors[side]) # Interpolates scale bar at the inttensity values for each fascicle in the paper
    FracList[np.where(FracList<0)] = 0

    hist = rv_histogram(np.histogram(FracList))
    frac = hist.rvs(size=1,random_state=rng)

    return frac*.01 # Converts from percentage to fraction

def getFiberTypeFractions(fascIdx, fascTypes):

    '''
    Finds the percentage of fibers of each type in the given fascicle, by sampling from the distribution in Figure 3 in Jayaprakash et al.
    '''

    rng = np.random.default_rng(fascIdx) # Sets random seed

    ###### In this block, we take the color intensity values from the plot of fiber type concentration per fascicle, and the scale bars.

    maffColors = [[103,211,191,255,254,191,157,254,231,232,255,248,255,226,244,255,227,212,192, 255],[160,80,82,81,83,82,118,158,135,182,184,182,134,110,102,128,82,107,119,114,135,126]] # Color intensities for each fascicle. First array corresponds to the left half of the nerve, second array to the right half of the nerve

    maffColorY = [0,10,15,20] # Fiber composition percentage from scale bar
    maffColorX = [255,205,151,95] # Corresponding color intensities for these composition percentages

    meffColors = [[243,253,180,169,226,202,169,214,198,207,219,177,180,210,226,171,169,169,169],[254,254,254,253,251,254,252,252,254,255,254,254,243,255,255,254,255,255,255,252,254]]

    meffColorY = [0,5,10,15]  # Fiber composition percentage from scale bar
    meffColorX = [255,233,208,185] # Corresponding color intensities for these composition percentages

    ueffColors = [[219,252,246,249,255,237,239,255,254,254,255,254,254,254,255,252,255,252,251,250,254],[196,210,210,197,220,195,211,227,232,241,246,233,248,242,234,220,222,235,236,193,194]]

    ueffColorY = [0,10,20,30]  # Fiber composition percentage from scale bar
    ueffColorX = [255,244,230,213]# Corresponding color intensities for these composition percentages

    ###########################

    if fascTypes[fascIdx]:
        side = 0
    else:
        side = 1

    if distribution_params['maff']['fiberTypeFractions'] is None:
        maffFrac = sampleFractionHistogram(maffColorX, maffColorY, maffColors, side, rng)
    else:
        maffFrac = distribution_params['maff']['fiberTypeFractions'][fascIdx] * 0.01

    if distribution_params['meff']['fiberTypeFractions'] is None:
        meffFrac = sampleFractionHistogram(meffColorX, meffColorY, meffColors, side, rng)
    else:
        meffFrac = distribution_params['meff']['fiberTypeFractions'][fascIdx] * 0.01

    ueffFrac = sampleFractionHistogram(ueffColorX,ueffColorY,ueffColors,side,rng)

    uaffFrac = 1-(maffFrac+meffFrac+ueffFrac)

    return maffFrac, meffFrac, ueffFrac, uaffFrac

def gammaDist(x,k,theta):

    return 1 / (gamma(k)*theta**k) * x**(k-1)*np.exp(-x/theta)

def prob(d, fiberType,diameter_params=None):

    '''
    Given a vector of fiber diameters d, and a particular fiber type, interpolates the fiber diameter probability distribution from the Jayaprakash paper over the vector d.
    '''

    d = d / pq.m # Removes units for compatibility reasons

    binSizeSamples = np.diff(d)[0]

    empiricalDiams, empiricalProbs = loadDiameterDistribution(fiberType)

    binSizeData = np.diff(empiricalDiams)[0] # Taking the first element ignores sloppy digitization towards the far end

    binRatio = binSizeSamples/binSizeData

    interp = interp1d(empiricalDiams,empiricalProbs,bounds_error=False,fill_value='extrapolate')

    interpD = interp(d)

    interpD[np.where(interpD<0)]=0

    if diameter_params is None:

        params = curve_fit(gammaDist,d*1e6,interpD*10,p0=[9,0.5],bounds=(0,np.inf)) # Fits gamma distribution to digitized data

    else:
        params = [[diameter_params[0],diameter_params[1]]]

    interpD = gammaDist(d*1e6,params[0][0],params[0][1]) * 0.1


    return (interpD * binRatio)/np.sum((interpD * binRatio).magnitude)


def MaffProb(d, maffProb,distribution_params,fascIdx=None):
    
    if distribution_params['maff']['diameterParams'] is not None:
        if fascIdx is not None and len(distribution_params['maff']['diameterParams'])>1:
            distributionParams = distribution_params['maff']['diameterParams'][fascIdx]
        else:
            distributionParams = distribution_params['maff']['diameterParams']
    else:
         distributionParams = None

    return maffProb * prob(d,'maff',distributionParams)

def MeffProb(d, meffProb,distribution_params,fascIdx=None):


    if distribution_params['meff']['diameterParams'] is not None:
        if fascIdx is not None and len(distribution_params['meff']['diameterParams'])>1:
            distributionParams = distribution_params['meff']['diameterParams'][fascIdx]
        else:
            distributionParams = distribution_params['meff']['diameterParams']
    else:
        distributionParams = None

    return meffProb * prob(d,'meff',distributionParams)

def UaffProb(d, uaffProb):

    uaffvals = np.loadtxt(r'D:\vagusOptimization\Data\uaffvals.csv',delimiter=',')

    return uaffProb * prob(d,'uaff')

def UeffProb(d, ueffProb):

    ueffvals = np.loadtxt(r'D:\vagusOptimization\Data\ueffvals.csv',delimiter=',')

    return ueffProb * prob(d,'ueff')

def getFasciclePositions():

    '''
    Loads spline positions from sim4life. For each fascicle, average spline positions to get fascicle position
    '''

    fasciclePositions = []

    positions = np.load(r'D:\vagusOptimization\Data\fiberPositions1950.npy',allow_pickle=True)

    pos = positions[0][1]

    for i in np.arange(1,len(positions)):
        pos = np.vstack((pos,positions[i][1]))

    for i in range(39): # Selects positions for each fascicle and averages them

        fiberPos = pos[i*50:(i+1)*50]

        fasciclePositions.append(np.mean(fiberPos,axis=0))

    return np.array(fasciclePositions)

def getFascicleTypes():

    fascPos = getFasciclePositions()

    # Selects whether fascicle should have more afferent or more efferent fibers, based on whether it is left or right (or above vs below) dividing line

    fascTypes = fascPos[:,0] > 8

    return fascTypes

def getFibersPerFascicle(fascIdx,fascTypes):

    '''
    Given the size of each fascicle, calculates the total number of fibers therein
    '''

    fascicleSizes = np.array([.24*.26,.16*.16,.18*.2,.16*.16,.12*.14,.16*.16,.1*.12,.24*.2,.2*.24,.18*.2,.14*.12,
                .16*.16,.1*.08,.16*.14,.12*.12,.08*.08,.14*.12,.1*.1,.2*.18,.14*.14,.14*.12,
                .12*.12,.22*.18,.14*.14,.14*.12,.18*.18,.16*.16,.1*.16,.12*.12,.22*.22,.1*.1,.1*.08,
                .12*.12,.1*.1,.12*.1,.14*.1,.1*.1,.14*.12,.18*.16])*1e-3**2


    fibersPerFascicle = getNumFibers(fascicleSizes,fascIdx,fascTypes)

    return fibersPerFascicle # Average value
