import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

from scipy.stats import rv_histogram

from math import gamma

from scipy.optimize import curve_fit

import quantities as pq

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.phiShape import *


def getAreaScaleFactor():
    
    ### Overall fiber counts in the nerve
    maffcount = 34000
    meffcount = 14800
    ueffcount = 21800 
    uaffcount = 315000
    ####
    
    #### Loads diameter distributions
    maffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/maffvals.csv',delimiter=',')
    meffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/meffvalsSmooth.csv',delimiter=',')
    uaffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/uaffvals.csv',delimiter=',')
    ueffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/ueffvals.csv',delimiter=',')
    #####
    
    #### Gets midpoints of histogram bins
    maffD = (maffvals[:-1,0] + maffvals[1:,0])/2 * 1e-6
    maffP = (maffvals[:-1,1] + maffvals[1:,1])/2 

    meffD = (meffvals[:-1,0] + meffvals[1:,0])/2 * 1e-6
    meffP = (meffvals[:-1,1] + meffvals[1:,1])/2

    uaffD = (uaffvals[:-1,0] + uaffvals[1:,0])/2 * 1e-6
    uaffP = (uaffvals[:-1,1] + uaffvals[1:,1])/2

    ueffD = (ueffvals[:-1,0] + ueffvals[1:,0])/2 * 1e-6
    ueffP = (ueffvals[:-1,1] + ueffvals[1:,1])/2
    #######
    
    
    maffArea = np.sum(maffD**2*maffP*maffcount / 100)
    meffArea = np.sum(meffD**2*meffP*meffcount / 100)
    uaffArea = np.sum(uaffD**2*uaffP*uaffcount / 100)
    ueffArea = np.sum(ueffD**2*ueffP*ueffcount / 100)

    totalFiberArea = maffArea + meffArea + uaffArea + ueffArea

    fascicleSizes = np.array([.24*.26,.16*.16,.18*.2,.16*.16,.12*.14,.16*.16,.1*.12,.24*.2,.2*.24,.18*.2,.14*.12,
                    .16*.16,.1*.08,.16*.14,.12*.12,.08*.08,.14*.12,.1*.1,.2*.18,.14*.14,.14*.12,
                .12*.12,.22*.18,.14*.14,.14*.12,.18*.18,.16*.16,.1*.16,.12*.12,.22*.22,.1*.1,.1*.08,
                .12*.12,.1*.1,.12*.1,.14*.1,.1*.1,.14*.12,.18*.16])*1e-3**2
    
    diamScaleFactor = (np.sum(fascicleSizes)/totalFiberArea)
    
    return diamScaleFactor

def getNumFibers(fascicleArea,diamScaleFactor,fascIdx,fascTypes):

    
    #### Assigns fiber counts per fascicle for each type
    
    maffFrac, meffFrac, ueffFrac, uaffFrac = getFiberTypeFractions(fascIdx,fascTypes)
        
    ###########
    
    ### Overall counts in the nerve
    maffcount = 34000
    meffcount = 14800
    ueffcount = 21800 
    uaffcount = 315000
    #######
    
    
    ##### Loads diameter distribution histograms
    maffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/maffvals.csv',delimiter=',')
    meffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/meffvalsSmooth.csv',delimiter=',')

    uaffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/uaffvals.csv',delimiter=',')
    ueffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/ueffvals.csv',delimiter=',')
    #######
    
    ### Gets midpoints of histogram bins
    maffD = (maffvals[:-1,0] + maffvals[1:,0])/2 *1e-6
    maffP = (maffvals[:-1,1] + maffvals[1:,1])/2 

    meffD = (meffvals[:-1,0] + meffvals[1:,0])/2 * 1e-6
    meffP = (meffvals[:-1,1] + meffvals[1:,1])/2

    uaffD = (uaffvals[:-1,0] + uaffvals[1:,0])/2 * 1e-6
    uaffP = (uaffvals[:-1,1] + uaffvals[1:,1])/2

    ueffD = (ueffvals[:-1,0] + ueffvals[1:,0])/2 * 1e-6
    ueffP = (ueffvals[:-1,1] + ueffvals[1:,1])/2
    #########

 
    maffArea = maffFrac * np.sum(maffD**2*maffP / 100)
    meffArea = meffFrac * np.sum(meffD**2*meffP / 100)
    uaffArea = uaffFrac * np.sum(uaffD**2*uaffP / 100)
    ueffArea = ueffFrac * np.sum(ueffD**2*ueffP / 100)

    fascicleNumber = fascicleArea[fascIdx] / (diamScaleFactor * (maffArea + meffArea + uaffArea + ueffArea)) 
    
#    fascicleNumber = 10000
    
    return maffFrac, meffFrac, uaffFrac, ueffFrac, fascicleNumber


def sampleFractionHistogram(ColorX,ColorY,Colors,side,rng):
    
    ##### In this function, we generate a histogram of per-fascicle fiber type fractions from the paper, and sample from it to produce the fraction for the simulated fascicle  
    
    Interp = interp1d(ColorX,ColorY,bounds_error=False,fill_value = 'extrapolate') # Interpolation object for the color scale bar
    FracList = Interp(Colors[side]) # Interpolates scale bar at the inttensity values for each fascicle in the paper
    FracList[np.where(FracList<0)] = 0
    
    hist = rv_histogram(np.histogram(FracList))
    frac = hist.rvs(size=1,random_state=rng)
    
    return frac*.01 # Converts from percentage to fraction

def getFiberTypeFractions(fascIdx, fascTypes):
    
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
        
    maffFrac = sampleFractionHistogram(maffColorX,maffColorY,maffColors,side,rng)
    meffFrac = sampleFractionHistogram(meffColorX,meffColorY,meffColors,side,rng)
    ueffFrac = sampleFractionHistogram(ueffColorX,ueffColorY,ueffColors,side,rng)
    
    uaffFrac = 1-(maffFrac+meffFrac+ueffFrac)
    
    return maffFrac, meffFrac, ueffFrac, uaffFrac

def gammaDist(x,k,theta):
    
    return 1 / (gamma(k)*theta**k) * x**(k-1)*np.exp(-x/theta)

def prob(d, vals,smooth):

    d = d / pq.m # Removes units for compatibility reasons
    
    binSizeSamples = np.diff(d)[0]
    
    empiricalDiams = vals[:,0]*1e-6 # From um to m
    empiricalProbs = vals[:,1]*0.01 # From percentage to fraction
    
    binSizeData = np.diff(empiricalDiams)[0] # Taking the first element ignores sloppy digitization towards the far end
    
    binRatio = binSizeSamples/binSizeData
    
    interp = interp1d(empiricalDiams,empiricalProbs,bounds_error=False,fill_value='extrapolate')
    
    interpD = interp(d)
    
    interpD[np.where(interpD<0)]=0
    
    if smooth:
        
        params = curve_fit(gammaDist,d*1e6,interpD*10,p0=[9,0.5])
        
        interpD = gammaDist(d*1e6,params[0][0],params[0][1]) * 0.1

                      
    return interpD * binRatio


def MaffProb(d, maffProb):
    
    maffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/maffvals.csv',delimiter=',')
    
    return maffProb * prob(d,maffvals,True)

def MeffProb(d, meffProb):
    
    meffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/meffvalsSmooth.csv',delimiter=',')
    
    return meffProb * prob(d,meffvals,True)

def UaffProb(d, uaffProb):
    
    uaffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/uaffvals.csv',delimiter=',')
    
    return uaffProb * prob(d,uaffvals,False)

def UeffProb(d, ueffProb):
    
    ueffvals = np.loadtxt('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/ueffvals.csv',delimiter=',')
    
    return ueffProb * prob(d,ueffvals,True)

def getFasciclePositions():
    
    '''
    Loads spline positions from sim4life. For each fascicle, average spline positions to get fascicle position
    '''
    
    fasciclePositions = []
    
    positions = np.load('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/fiberPositions1950.npy',allow_pickle=True)
    
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
    
    diamScaleFactor = getAreaScaleFactor()
        
    fascicleSizes = np.array([.24*.26,.16*.16,.18*.2,.16*.16,.12*.14,.16*.16,.1*.12,.24*.2,.2*.24,.18*.2,.14*.12,
                .16*.16,.1*.08,.16*.14,.12*.12,.08*.08,.14*.12,.1*.1,.2*.18,.14*.14,.14*.12,
                .12*.12,.22*.18,.14*.14,.14*.12,.18*.18,.16*.16,.1*.16,.12*.12,.22*.22,.1*.1,.1*.08,
                .12*.12,.1*.1,.12*.1,.14*.1,.1*.1,.14*.12,.18*.16])*1e-3**2
    
    maffFrac, meffFrac, uaffFrac, ueffFrac, fibersPerFascicle = getNumFibers(fascicleSizes,diamScaleFactor,fascIdx,fascTypes)
    
    return fibersPerFascicle # Average value