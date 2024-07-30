import numpy as np
import pandas as pd

from mpi4py import MPI

import quantities as pq

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *

def getFascIdx():

    ## Selects fascicle and random seed for each available cpu
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()

    if comm.Get_size() != 39:
        raise AssertionError('Must use exactly 39 ranks')
    
    fascIdx = int(rank % 39)
    #####

    return fascIdx

def loadActionPotentialShapes():

    ### Loads action potential shapes in time
    ap = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShape20.xlsx') # Rat
    ap2 = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapePoint8.xlsx') # Sundt
    ####

    aps = [ap,ap2]

    return aps

def getTime():

    nx=500000

    tmin=-3 # In s
    tmax=3 # In s
    time=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))*pq.s

    return time

def convolveToGetSignal(time, current, phi, recordingCurrent, variance=np.array([0])):

    aps = loadActionPotentialShapes()

    Vs = getVs(aps,time) # Interpolates action potential shape in time

    signals = []

    for i, V in enumerate(Vs):  

        der = np.diff(V,n=2)/((time[1]-time[0])**2) # Second derivative of action potential shape

        cv = []
        
        for j in range(len(current)):

            for k in range(len(variance)):
            
                c = fftconvolve(der,phi[i,:,j,k],mode='same') # Convolves second derivative with exposure
                
                cv.append(c)

        cv = np.array(cv)

        signals.append(cv)

    signals = np.array(signals)
    
    signals /= recordingCurrent

    return signals

def saveSignals(outputfolder, distanceIdx, fascIdx, signals):

    names = ['maff','meff','uaff','ueff']

    for typeIdx in range(len(names)):

        np.save(outputfolder+'/'+names[typeIdx]+'/'+str(distanceIdx)+'/signals_'+str(fascIdx)+'.npy',signals[typeIdx])

def getDiameterScalingOfCurrent(d, time, velocityList):

    '''
    For each diameter, returns scaling factor eta, which links the transmembrane current to the diameter of the fiber
    '''

    scaling0 = Scaling(d,0)* (time[1]-time[0])/velocityList[0]
    scaling1 = Scaling(d,1)* (time[1]-time[0])/velocityList[1]

    return scaling0, scaling1

def getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, outputfolder, variance=np.array([0])):

    phiWeightMaff, phiWeightMeff, phiWeightUaff, phiWeightUeff = getPhiWeight(d,current,fascIdx, fascTypes, stimulusDirectory, variance) # For each of the four fiber types, returns scaling factor for each diameter

    myelinatedCurrentScaling, unmyelinatedCurrentScaling = getDiameterScalingOfCurrent(d, time, velocityList)

    myelinatedCurrentScaling = myelinatedCurrentScaling.magnitude
    unmyelinatedCurrentScaling = unmyelinatedCurrentScaling.magnitude

    phiWeightMaff = phiWeightMaff.magnitude
    phiWeightMeff = phiWeightMeff.magnitude
    phiWeightUaff = phiWeightUaff.magnitude
    phiWeightUeff = phiWeightUeff.magnitude
    
   #np.save(outputfolder+'/diameters/scaling_'+str(fascIdx)+'.npy',[scaling0,scaling1])
    
    maffscaling = phiWeightMaff.T * myelinatedCurrentScaling[:,np.newaxis] 
    meffscaling = phiWeightMeff.T * myelinatedCurrentScaling[:,np.newaxis]
    uaffscaling = phiWeightUaff.T * unmyelinatedCurrentScaling[:,np.newaxis]
    ueffscaling = phiWeightUeff.T * unmyelinatedCurrentScaling[:,np.newaxis]
    
   # np.save(outputfolder+'/maff/scaling_'+str(fascIdx)+'.npy',maffscaling)
    
    #np.save(outputfolder+'/meff/scaling_'+str(fascIdx)+'.npy',meffscaling)
    
   # np.save(outputfolder+'/uaff/scaling_'+str(fascIdx)+'.npy',uaffscaling)
    #
   # np.save(outputfolder+'/ueff/scaling_'+str(fascIdx)+'.npy',ueffscaling)

    return [maffscaling, meffscaling, uaffscaling, ueffscaling]

def getExposureFunctions(phiShapesByType, scalingFactorsByType, outputfolder, distanceIdx):

    phi = [0,0,0,0]

    
    phiShapeMyelinated, phiShapeUnmyelinated = phiShapesByType

    maffScaling, meffScaling, uaffScaling, ueffScaling = scalingFactorsByType

    
    phi[0] += phiShapeMyelinated.T @ maffScaling 
    
    phi[1] += phiShapeMyelinated.T @ meffScaling 
    
    phi[2] += phiShapeUnmyelinated.T @ uaffScaling 
    
    phi[3] += phiShapeUnmyelinated.T @ ueffScaling 

    phi = np.array(phi)
    
   # np.save(outputfolder+'/phis/'+str(distanceIdx)+'/'+str(fascIdx)+'.npy',phi)

    return phi

def getPhiShapes(fascIdx, distance, recordingDirectory, velocityList, time):
    
    phiFunc = FitPhiShape(fascIdx, distance, recordingDirectory)# Defines an interpolation function for the recording exposure for the fasicle
    
    ### For each diameter, defines a shifted and scaled exposure function
    phiShapeMyelinated = PhiShape(velocityList[0],time,phiFunc)
    
    phiShapeUnmyelinated = PhiShape(velocityList[1],time,phiFunc)

    #np.save(outputfolder+'/phis/'+str(distanceIdx)+'/rawShapes_'+str(fascIdx)+'.npy',[phiShape0,phiShape1])

    return [phiShapeMyelinated, phiShapeUnmyelinated]

def getDistance(distanceIdx, recording):

    if 'distances' in recording.keys():
        distances = recording['distances']*pq.m
    else:
        distances = [0.06,0.01]*pq.m # Stimulus-recording distance, in m
    
    distance = distances[distanceIdx]

    return distance

def getVariance(stimulus):

    if 'variance' in stimulus.keys():

        variance = stimulus['variance']
    else:
        variance=np.array([0])
    
def runSim(outputfolder, distanceIdx, stimulus, recording):
   
    fascIdx = getFascIdx()
    
    current = stimulus['current'] # Current applied in finite element simulation of recruitment
    stimulusDirectory = stimulus['stimulusDirectory'] # Location of titration outputs from S4:
    
    distance = getDistance(distanceIdx, recording)

    time = getTime()

    recordingCurrent = recording['recordingCurrent']*pq.A # Current in the S4L recording simulation
    recordingDirectory = recording['recordingDirectory'] # Location of exported potential fields interpolated along fascicle centers

    d = getDiameters() 

    velocityList = getVelocities(d) # Gets velocity for each diameter
    
    fascTypes = getFascicleTypes() # Defines whether fasicle is on left or right side of nerve  

    scalingFactorsByType = getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, outputfolder, variance)

    phiShapesByType = getPhiShapes(fascIdx, distance, recordingDirectory, velocityList, time)
        
    phi = getExposureFunctions(phiShapesByType, scalingFactorsByType, outputfolder, distanceIdx)
    
    signals = convolveToGetSignal(time, current, phi, recordingCurrent, variance)

    saveSignals(outputfolder, distanceIdx, fascIdx, signals)
    
