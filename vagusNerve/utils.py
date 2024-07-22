import numpy as np
import pandas as pd


from scipy.stats import norm
from scipy.io import loadmat

from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.stats import norm
from scipy.fft import fft, ifft, fftshift,ifftshift
from scipy.signal import fftconvolve, butter, sosfilt

from scipy.stats import rv_histogram
from mpi4py import MPI

from math import gamma

from scipy.optimize import curve_fit

import quantities as pq

import sys

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *

def Scaling(d,fiberType): # Diameter dependent scaling of signal
    
    if fiberType == 0: # Myelinated fiber
        resistivity_intracellular = 0.7  # ohm meters
        deff = d # Calculates internal diameter from external diameter
        
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



def getVs(aps,tphi): # Interpolates AP shapes in time for myelinated and unmyelinated fibers
    
    Vs = []

    for i, ap in enumerate(aps): # Iterates throguh the two classes (meylinated and unlyeminated)

        v = FitAPShape(ap,tphi)

        Vs.append(v)
        Vs.append(v)
        
    return Vs