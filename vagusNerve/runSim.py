# Apache-2.0

import numpy as np
import pandas as pd


import quantities as pq

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *

import time as tm

def loadActionPotentialShapes():

    ### Loads action potential shapes in time
    ap = pd.read_excel(r'D:\vagusOptimization\Data\APShape20.xlsx') # Rat
    ap2 = pd.read_excel(r'D:\vagusOptimization\Data\APShapePoint8.xlsx') # Sundt
    ####

    aps = [ap,ap2]

    return aps

def getTime():

    nx=50000

    tmin=-.5 # In s
    tmax=.5 # In s
    time=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))*pq.s

    return time

def convolveToGetSignal(time, current, phi, recordingCurrent, variance=np.array([0])):

    aps = loadActionPotentialShapes()

    Vs = getVs(aps,time) # Interpolates action potential shape in time

    V = Vs[0]

    signals = []

    der = np.diff(V,n=2)/((time[1]-time[0])**2) # Second derivative of action potential shape

    cv = []

    for i in range(1):

        for j in range(np.max( (len(current),len(variance)) )):


            c = fftconvolve(der,phi[i,:,j],mode='same') # Convolves second derivative with exposure

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

    return scaling0

def getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, distribution_params, variance=np.array([0])):

    phiWeightMaff, phiWeightMeff = getPhiWeight(d,current,fascIdx, fascTypes, stimulusDirectory, distribution_params, variance) # For each of the four fiber types, returns scaling factor for each diameter

    myelinatedCurrentScaling = getDiameterScalingOfCurrent(d, time, velocityList)

    myelinatedCurrentScaling = myelinatedCurrentScaling.magnitude

    phiWeightMaff = phiWeightMaff.magnitude
    phiWeightMeff = phiWeightMeff.magnitude


    maffscaling = phiWeightMaff.T * myelinatedCurrentScaling[:,np.newaxis]
    meffscaling = phiWeightMeff.T * myelinatedCurrentScaling[:,np.newaxis]


    return [maffscaling, meffscaling]

def getExposureFunctions(phiShapesByType, scalingFactorsByType, distanceIdx, fascIdx):



    phiShapeMyelinated = phiShapesByType

    maffScaling, meffScaling = scalingFactorsByType


    phi0 = phiShapeMyelinated.T @ maffScaling

    phi1 = np.zeros_like(phi0)#phiShapeMyelinated.T @ meffScaling


    phi = np.array([phi0,phi1])

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

def runSim(distanceIdx, stimulus, recording, fascIdx, distribution_params, numDiameters=2000):

    t = tm.time()
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
    
    scalingFactorsByType = getScalingFactors(d,current,fascIdx, fascTypes, stimulusDirectory, time, velocityList, distribution_params, variance)
    
    cutoff = getPhiCutoff(recording)

    phiShapesByType = getPhiShapes(fascIdx, distance, recordingDirectory, velocityList, time, cutoff)

    phi = getExposureFunctions(phiShapesByType, scalingFactorsByType, distanceIdx, fascIdx)
    signals = convolveToGetSignal(time, current, phi, recordingCurrent, variance)
    return signals
