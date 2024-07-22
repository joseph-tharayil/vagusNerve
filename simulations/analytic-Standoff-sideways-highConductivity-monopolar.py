import numpy as np
import pandas as pd


from scipy.stats import norm
from scipy.io import loadmat

from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.stats import norm
import multiprocessing as mp
from scipy.fft import fft, ifft, fftshift,ifftshift
from scipy.signal import fftconvolve, butter, sosfilt

from scipy.stats import rv_histogram
from mpi4py import MPI

from math import gamma

from scipy.optimize import curve_fit

import quantities as pq

import sys

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
    
#         N = 5
#         empiricalDiams = np.convolve(empiricalDiams, np.ones(N)/N, mode='valid')
#         empiricalProbs = np.convolve(empiricalProbs, np.ones(N)/N, mode='valid')
    
                      
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


def sortTitrationSpace(table): 
    
    '''
    This function loads the titration results table from Sim4Life, and sorts it such that the splines and fascicles are in numerical order
    '''
    
    fascicles = []
    splines = []
    
    for column in table.columns:
                
        fasc = column.split('_')[0]
        
        try: # Get fascicle and spline name from column name
            fasc = int(fasc.split(' ')[-1])
            fascicles.append(fasc)
            spline = column.split(' [')[0].split('_')[-1]
            splines.append(int(spline))
        except: # Fasicle 0 is just called Fascicle_, not Fascicle_i like the other ones
            fascicles.append(0)
            spline = column.split(' [')[0].split('_')[-1]
            splines.append(int(spline))

    fascicles = np.array(fascicles)
    splines = np.array(splines)

    indices = np.lexsort((splines,fascicles))
    
    return table.iloc[:,indices]

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
    
def getFascicleTypes(iteration):
    
    fascPos = getFasciclePositions()
    
    nerveCenter = np.mean(fascPos,axis=0)
    
    # Selects whether fascicle should have more afferent or more efferent fibers, based on whether it is left or right (or above vs below) dividing line
    
    if iteration == 0:
        fascTypes = fascPos[:,0] > 8
    elif iteration == 1:
        fascTypes = fascPos[:,0] < 8
    elif iteration == 2:
        fascTypes = fascPos[:,1] > -9
    else:
        fascTypes = fascPos[:,1] < -9
    
    return fascTypes

def getFibersPerFascicle(fascIdx,fascTypes):
    
    diamScaleFactor = getAreaScaleFactor()
        
    fascicleSizes = np.array([.24*.26,.16*.16,.18*.2,.16*.16,.12*.14,.16*.16,.1*.12,.24*.2,.2*.24,.18*.2,.14*.12,
                .16*.16,.1*.08,.16*.14,.12*.12,.08*.08,.14*.12,.1*.1,.2*.18,.14*.14,.14*.12,
                .12*.12,.22*.18,.14*.14,.14*.12,.18*.18,.16*.16,.1*.16,.12*.12,.22*.22,.1*.1,.1*.08,
                .12*.12,.1*.1,.12*.1,.14*.1,.1*.1,.14*.12,.18*.16])*1e-3**2
    
    maffFrac, meffFrac, uaffFrac, ueffFrac, fibersPerFascicle = getNumFibers(fascicleSizes,diamScaleFactor,fascIdx,fascTypes)
    
    return fibersPerFascicle # Average value

def Recruitment(current,diameters, fascIdx):
    
    d0Myelinated = 4e-6
    d0Unmyelinated = 0.8e-6
    
   
    #### Loads and sorts titration factors from S4L. Sorting is justg to make sure that fibers and fascicles are in numerical order (ie, fiber 0-fiber50, fascicle0-fascicle39)
    
    titrationFactorsMeff = sortTitrationSpace(pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/Titration_SENN_Sideways.xlsx',index_col=0)).iloc[-1].values
    
    titrationFactorsUaff = sortTitrationSpace(pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_HighConductivity.xlsx',index_col=0)).iloc[-1].values


    ####


    titrationFactors = [titrationFactorsMeff, titrationFactorsUaff]


    for j in [0,1]: # Myelinated and unmyelinated, respectively

        titrationFac = np.array(titrationFactors[j][fascIdx*50:(fascIdx+1)*50]) # Selects fibers in fascicle


        if j == 0:

            myelinatedTitrationFactors = titrationFac

        else:

            unmyelinatedTitrationFactors = titrationFac

### Defines CDF of the titration curves
#     hist, edges = np.histogram(myelinatedTitrationFactors,bins=10,density=True)
#     cdf = np.cumsum(hist)/np.sum(hist)
#     midpts = (edges[:-1]+edges[1:])/2

    cdf = np.arange(len(myelinatedTitrationFactors))/len(myelinatedTitrationFactors)
    midpts = np.sort(myelinatedTitrationFactors)
    dupIdx =  np.where(np.diff(midpts)==0)
    midpts = np.delete(midpts,dupIdx)
    cdf = np.delete(cdf,dupIdx)
    
    interp = interp1d(midpts,cdf,bounds_error=False,fill_value=(0,1))
    
#     histU, edgesU = np.histogram(unmyelinatedTitrationFactors,bins=10,density=True)
#     cdfU = np.cumsum(histU)/np.sum(histU)
#     midptsU = (edgesU[:-1]+edgesU[1:])/2
#     interpU = interp1d(midptsU,cdfU,bounds_error=False,fill_value=(0,1))

    cdfU = np.arange(len(unmyelinatedTitrationFactors))/len(unmyelinatedTitrationFactors)
    midptsU = np.sort(unmyelinatedTitrationFactors)
    dupIdxU =  np.where(np.diff(midptsU)==0)
    midptsU = np.delete(midptsU,dupIdxU)
    cdfU = np.delete(cdfU,dupIdxU)
    
    interpU = interp1d(midptsU,cdfU,bounds_error=False,fill_value=(0,1))
##############
        
    myelinated = []
    unmyelinated = []


    for j, d in enumerate(diameters):

        myelinated.append(interp(current*(d/d0Myelinated)))

        unmyelinated.append(interpU(current*d/d0Unmyelinated))

    
    return [myelinated,unmyelinated]

def Scaling(d,fiberType): # Diameter dependent scaling of signal
    
    if fiberType == 0: # Myelinated fiber
        resistivity_intracellular = 1.1  # ohm meters
        deff = 0.7 * d# 0.32 * d # Calculates internal diameter from external diameter
        
    elif fiberType == 1: # Unmyelinated Fiber, Sundt Model
        resistivity_intracellular = 1
        deff = d
        
    elif fiberType == 2: # Unmyelinated fiber, tigerholm model
        resistivity_intracellular = 0.354 # ohm meters
        deff = d
    
    segment_length = 50e-6
    
    surfaceArea = segment_length * deff 
    
    xSectionArea = np.pi * (deff/2)**2
        
    resistance_intracellular = resistivity_intracellular/xSectionArea
        
    
    current_scale_factor = 1/(resistance_intracellular) 
                         
             
    return current_scale_factor


def editPhiShape(phi,distance):
    
    ''' 
    This function takes the recording exposure curve from S4L, shifts it to match the desired distance from stimulus to recording, and smooths it
    '''
    
    xvals = phi.iloc[:,0].values+distance -phi.iloc[np.argmax(phi.iloc[:,1].values),0] # Shift to match desired distance

    phiShapeEmpirical = phi.iloc[:,1].values-np.mean(phi.iloc[:,1])

    
   ######## 
    
    ####### Makes sure that the potential at the proximal end of the fiber goes all the way to zero

    if np.any(phiShapeEmpirical[:np.argmax(phiShapeEmpirical)]<0): # If the potential is negative at the end of the fiber, sets the potential to 0
        
        
        first = np.where(phiShapeEmpirical[:np.argmax(phiShapeEmpirical)]<0)[0][-1]    
        
        phiShapeEmpirical[:first] = 0
        
    else: # If the potential does not go all the way to 0 by the end of the fiber, forces it to zero
        
        first = np.where(np.abs(np.diff(phiShapeEmpirical))>5e-6)[0][0] # Based on derivative of function, selects point after whcih not to change values
        

        ### Linearizes potential before this point, up until it reaches 0 
        firsta = np.where(phiShapeEmpirical[first]-5e-5*np.arange(first)<0)[0][0]
        firsta = first-firsta

        phiShapeEmpirical[firsta:first] = 5e-5*np.arange(first-firsta)
        #######
        
        phiShapeEmpirical[0:firsta]=0 # Sets potential to zero
        
    ############
   
    #### Does the same kind of smoothing as above, but for the distal end fo the fiber
    if np.any(phiShapeEmpirical[np.argmax(phiShapeEmpirical):]<0):
                
        last = np.where(phiShapeEmpirical[np.argmax(phiShapeEmpirical):]<0)[0][0]+np.argmax(phiShapeEmpirical)
        
        phiShapeEmpirical[last:] = 0
        
#     else:
#         last = np.where(np.abs(np.diff(phiShapeEmpirical))>5e-5)[0][-1]
#         lasta = np.where(phiShapeEmpirical[last]+5e-6*np.arange(len(phiShapeEmpirical)-last)>0)[0][0]
#         lasta += last

#         phiShapeEmpirical[last:lasta] = 5e-5* np.arange(lasta-last)+ phiShapeEmpirical[last]
#         phiShapeEmpirical[lasta:] = 0
    

    return xvals, phiShapeEmpirical

def FitPhiShape(fascIdx,distance):
    
    ''' 
    This function creates an interpolation object for the recording exposure
    '''

    phi = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/PhiConductivity_Monopolar_Corrected/'+str(fascIdx)+'_BetterConductivity.xlsx')
    
    phi = phi.iloc[:14000]
    
    xvals, phiShapeEmpirical = editPhiShape(phi,distance)

    return interp1d(xvals,phiShapeEmpirical,bounds_error=False,fill_value=(phiShapeEmpirical[0],phiShapeEmpirical[-1]))

def PhiShape(velocity,t,function):
    
    
    '''
    This function stretches the recording exposure in time, based on fiber velocity
    '''

    phiOut = []
    
#     sos = butter(1, 20000, 'lp', fs=83333, output='sos') # We apply a low-pass filter to avoid numerical issues

    for i, v in enumerate(velocity): # Iterates through fibers

        x = t*v ### Defines interpolation points from time vector and fiber velocity

        out = function(x)
            
#         filtered = sosfilt(sos,out)


        phiOut.append(out)

    return np.array(phiOut)

def getVelocities(d0List,velocities,dList):

    velocityList = []

    for i, d in enumerate(d0List):
        
        if i < 1: # Myelinated velocities are linear with diameter
            d0 = d0List[0]
            v0 = velocities[0]
            
            velocity = v0 * dList/d0
            
        else: # Unmyelinated velocities go as the square root of diameter
            d0 = d0List[1]
            v0 = velocities[1]
            
            velocity = v0 * np.sqrt(dList/d0)

        
        velocityList.append(velocity)

    return velocityList


def PhiWeight(d, current,fascIdx, fascTypes):
    
    phiWeight = [ [[],[]], [[],[]] ]
    
    recruitment = Recruitment(current,d,fascIdx)
    
    scaling = []
    
    scalingFactors = [1,2]
    
    maffProb, meffProb, ueffProb, uaffProb = getFiberTypeFractions(fascIdx, fascTypes)
    
    numFibersPerFascicle = getFibersPerFascicle(fascIdx,fascTypes)
    
    
##### Weight is given by the product of the recruitment curve and the diameter probability curve
    phiWeight[0][0] =  MaffProb(d,maffProb)  * recruitment[0] * numFibersPerFascicle
    phiWeight[0][1] =  MeffProb(d,meffProb)  * recruitment[0] * numFibersPerFascicle
    
    
    phiWeight[1][0] =  UaffProb(d,uaffProb)  * recruitment[-1] * numFibersPerFascicle
    phiWeight[1][1] =  UeffProb(d,ueffProb)  * recruitment[-1] * numFibersPerFascicle
    
    
    np.save(outputfolder+'/'+str(0)+'/fascicles'+'/fibers'+str(fascIdx)+'.npy',numFibersPerFascicle)
    np.save(outputfolder+'/'+str(0)+'/fascicles'+'/probs'+str(fascIdx)+'.npy',[maffProb,meffProb,uaffProb,ueffProb])
    np.save(outputfolder+'/'+str(0)+'/fascicles'+'/probDist'+str(fascIdx)+'.npy',[MaffProb(d,maffProb),MeffProb(d,meffProb),UaffProb(d,maffProb),UeffProb(d,meffProb)])
        
    
    return phiWeight,recruitment


def FitAPShape(ap,tphi): # Interpolates AP shape for a given AP
    
    
    # Ignores initial transient
    tv = ap.iloc[50:,0]
    v = ap.iloc[50:,1]
    
    ### Sets peak time to 0
    peak = tv[np.argmax(v)]
    tv -= peak


    apShapeEmpirical = v.values
        

    func = interp1d(tv,apShapeEmpirical,bounds_error=False,fill_value=(apShapeEmpirical[0],apShapeEmpirical[-1]))
    
    Vs = func(tphi)  
    
    
    #### Applies low-pass filter with very high cutoff, to remove artifacts
    sos = butter(1, 20000, 'lp', fs=83333, output='sos')
    
    V = sosfilt(sos,Vs)
    
    V[:10] = V[10]
         
    return V

def getDiameters(iteration):
    
   
    minDiam = .1
    
    
    maxDiam = 15 #7 + 5*iteration/30 
    
    d = np.linspace(minDiam,maxDiam,2000)*1e-6

    return d

def getPhiWeight(d, current,fascIdx,fascTypes):
    
    phiWeight = []
    recruitment = []
    
    for c in current:
        
        p, rec = PhiWeight(d,c,fascIdx,fascTypes)

        phiWeight.append(p)
        recruitment.append(rec)
        
    np.save(outputfolder+'/'+str(0)+'/recruitment'+'/recruitment_'+str(fascIdx)+'.npy',recruitment)

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

def getVs(aps,tphi): # Interpolates AP shapes in time for myelinated and unmyelinated fibers
    
    Vs = []

    for i, ap in enumerate(aps): # Iterates throguh the two classes (meylinated and unlyeminated)

        v = FitAPShape(ap,tphi)

        Vs.append(v)
        Vs.append(v)
        
    return Vs

    
def main(outputfolder,distanceIdx):
    
   
    ## Selects fascicle and random seed for each available cpu
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    iteration = int(rank/39)
    fascIdx = int(rank % 39)
    #####
    
    distances = [0.06,.05,0.01] # Stimulus-recording distance, in m
    
    #### Random generators
    diameterGenerator = np.random.default_rng(seed=fascIdx)# Diameters depend on fascicle
    latencyGenerator = np.random.default_rng(seed=iteration) # Latency depends on trial index
    #####
       
    
    distance = distances[distanceIdx]
    
    ### Loads action potential shapes in time
    ap = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapeSenn.xlsx') # Rat
    ap2 = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapePoint8.xlsx') # Sundt
    ####

    aps = [ap,ap2]

    nx=500000

    tmin=-3 # In s
    tmax=3 # In s
    tphi=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))

    current = np.linspace(100,500,5)/28.6
#    current = np.linspace(100,500,5)/241 # Currents from 100 to 1000 uA, scaled by the current in the S4L titration simulation
    
#     current = np.linspace(100,500,5)/180 # Currents from 100 to 1000 uA, scaled by the current in the S4L titration simulation
    
    recordingCurrent = 460e-6 # Current in the S4L recording simulation

    d0List = [20*1e-6,.8*1e-6] # Diameters of myelinated and unmyelinated fibers used to calculate velocities

    velocities = [84.42,0.416] # Velcities for the above diamters
    
    
    fascTypes = getFascicleTypes(iteration) # Defines whether fasicle is on left or right side of nerve

    d = getDiameters(iteration)   


    phiWeight0, phiWeight1, phiWeight2, phiWeight3 = getPhiWeight(d,current,fascIdx, fascTypes) # Scaling factor for each diameter


    phi = [0,0,0,0]


    velocityList = getVelocities(d0List,velocities,d) # Gets velocity for each diameter


    names = ['myelinated','unmyelinated']


    phiFunc = FitPhiShape(fascIdx,distance)# Defines an interpolation function for the recording exposure for the fasicle
    
    ### For each diameter, defines a shifted and scaled exposure function
    phiShape0 = PhiShape(velocityList[0],tphi,phiFunc)
    
    phiShape1 = PhiShape(velocityList[1],tphi,phiFunc)
    
    np.save(outputfolder+'/'+str(iteration)+'/'+'phis/'+str(distanceIdx)+'/rawShapes_'+str(fascIdx)+'.npy',[phiShape0,phiShape1])
   #############


### Scales exposure functions

    scaling0 = (Scaling(d,0)* (tphi[1]-tphi[0])/velocityList[0])
    scaling1 = (Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])
    
    np.save(outputfolder+'/'+str(iteration)+'/diameters/scaling_'+str(fascIdx)+'.npy',[scaling0,scaling1])
    
    scaling00 = phiWeight0.T*scaling0[:,np.newaxis]
    scaling10 = phiWeight1.T*scaling0[:,np.newaxis]
    scaling01 = phiWeight2.T*scaling1[:,np.newaxis]
    scaling11 = phiWeight3.T*scaling1[:,np.newaxis]
    
    np.save(outputfolder+'/'+str(iteration)+'/maff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling00)
    
    np.save(outputfolder+'/'+str(iteration)+'/meff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling10)
    
    np.save(outputfolder+'/'+str(iteration)+'/uaff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling01)
    
    np.save(outputfolder+'/'+str(iteration)+'/ueff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling11)
    
    phi[0] += np.matmul(phiShape0.T,scaling00)
    
    phi[1] += np.matmul(phiShape0.T,scaling10)
    
    phi[2] += np.matmul(phiShape1.T,phiWeight2.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])
    
    phi[3] += np.matmul(phiShape1.T,phiWeight3.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])


    phi = np.array(phi)
    
    np.save(outputfolder+'/'+str(0)+'/'+'phis/'+str(distanceIdx)+'/'+str(fascIdx)+'.npy',phi)
################


    signal = 0

    Vs = getVs(aps,tphi) # Interpolates action potential shape in time

    k = np.linspace(-1,1-(2.0/np.shape(Vs)[1]),np.shape(Vs)[1])

    signals = []

    for i, V in enumerate(Vs):  

        der = np.diff(V,n=2)/((tphi[1]-tphi[0])**2) # Second derivative of action potential shape

        cv = []
        
        for j in range(len(current)):
            
            c = fftconvolve(der,phi[i,:,j],mode='same') # Convolves second derivative with exposure
            
            cv.append(c)

        cv = np.array(cv)

        signals.append(cv)

    signals = np.array(signals)
    
    signals /= recordingCurrent


    names = ['maff','meff','uaff','ueff']

    for typeIdx in range(len(names)):

        np.save(outputfolder+'/'+str(iteration)+'/'+names[typeIdx]+'/'+str(distanceIdx)+'/signals_'+str(fascIdx)+'.npy',signals[typeIdx])


if __name__=="__main__":
    
    
    outputfolder = sys.argv[1]
    distanceIdx = int(sys.argv[2])
    
    main(outputfolder,distanceIdx)
