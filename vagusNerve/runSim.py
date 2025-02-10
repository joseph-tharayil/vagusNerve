# Apache-2.0

import numpy as np
import pandas as pd


import quantities as pq

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *


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

    V = Vs[0]


    der = np.diff(V,n=2)/((time[1]-time[0])**2) # Second derivative of action potential shape

    cv = []

    for j in range(np.max( (len(current),len(variance)) )):
        cv = []

        c = fftconvolve(der,phi[:,j],mode='same') # Convolves second derivative with exposure

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

def getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, distribution_params, variance=np.array([0])):

    phiWeightMaff, phiWeightMeff, phiWeightUaff, phiWeightUeff = getPhiWeight(d,current,fascIdx, fascTypes, stimulusDirectory, distribution_params, variance) # For each of the four fiber types, returns scaling factor for each diameter

    myelinatedCurrentScaling, unmyelinatedCurrentScaling = getDiameterScalingOfCurrent(d, time, velocityList)

    myelinatedCurrentScaling = myelinatedCurrentScaling.magnitude
    unmyelinatedCurrentScaling = unmyelinatedCurrentScaling.magnitude

    phiWeightMaff = phiWeightMaff.magnitude
    phiWeightMeff = phiWeightMeff.magnitude
    phiWeightUaff = phiWeightUaff.magnitude
    phiWeightUeff = phiWeightUeff.magnitude

    maffscaling = phiWeightMaff.T * myelinatedCurrentScaling[:,np.newaxis]
    meffscaling = phiWeightMeff.T * myelinatedCurrentScaling[:,np.newaxis]
    uaffscaling = phiWeightUaff.T * unmyelinatedCurrentScaling[:,np.newaxis]
    ueffscaling = phiWeightUeff.T * unmyelinatedCurrentScaling[:,np.newaxis]

    return [maffscaling, meffscaling, uaffscaling, ueffscaling]

def getExposureFunctions(phiShapesByType, scalingFactorsByType, distanceIdx, fascIdx):

    phiShapeMyelinated = phiShapesByType

    maffScaling, meffScaling, uaffScaling, ueffScaling = scalingFactorsByType

    phi0 = phiShapeMyelinated.T @ maffScaling

    phi1 = phiShapeMyelinated.T @ meffScaling

    phi = phi0+phi1


    return phi

def getPhiShapes(fascIdx, distance, recordingDirectory, velocityList, time, cutoff=1e-4):

    phiFunc = FitPhiShape(fascIdx, distance, recordingDirectory, cutoff)# Defines an interpolation function for the recording exposure for the fasicle

    ### For each diameter, defines a shifted and scaled exposure function
    phiShapeMyelinated = PhiShape(velocityList[0],time,phiFunc)

    return phiShapeMyelinated

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
    return variance

def getPhiCutoff(recordingDirectory):

    '''
    Slope at which to linearize the exposure curve obtained from the FEM simulation. Defaults to 1e-4.
    Used by the FitPhiShape function
    '''

    if 'cutoff' in recordingDirectory.keys():

        cutoff = recordingDirectory['cutoff']
    else:
        cutoff = 1e-4

    return cutoff

def runSim(fascIdx,stimulus=None,recording=None,numDiameters=2000):

    distanceIdx = 0

    current = stimulus['current'] # Current applied in finite element simulation of recruitment
    stimulusDirectory = stimulus['stimulusDirectory'] # Location of titration outputs from S4:

    distance = getDistance(distanceIdx, recording)

    time = getTime()

    recordingCurrent = recording['recordingCurrent']*pq.A # Current in the S4L recording simulation
    recordingDirectory = recording['recordingDirectory'] # Location of exported potential fields interpolated along fascicle centers

    d = getDiameters(numDiameters)

    velocityList = getVelocities(d) # Gets velocity for each diameter

    fascTypes = getFascicleTypes() # Defines whether fasicle is on left or right side of nerve

    variance = getVariance(stimulus)

    if len(variance) > 1 and len(current)>1:
        raise AssertionError('Either variance or current must be constant')

    scalingFactorsByType = getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, outputfolder, variance)


    cutoff = getPhiCutoff(recording)

    phiShapesByType = getPhiShapes(fascIdx, distance, recordingDirectory, velocityList, time, cutoff)

    phi = getExposureFunctions(phiShapesByType, scalingFactorsByType, outputfolder, distanceIdx, fascIdx)


    signals = convolveToGetSignal(time, current, phi, recordingCurrent, variance)

    return signals
