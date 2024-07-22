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
        
        first = np.where(np.abs(np.diff(phiShapeEmpirical))>1e-4)[0][0] # Based on derivative of function, selects point after whcih not to change values
        

        ### Linearizes potential before this point, up until it reaches 0 
        firsta = np.where(phiShapeEmpirical[first]-1e-4*np.arange(first)<0)[0][0]
        firsta = first-firsta

        phiShapeEmpirical[firsta:first] = 1e-4*np.arange(first-firsta)
        #######
        
        phiShapeEmpirical[0:firsta]=0 # Sets potential to zero
        
    ############
   
    #### Does the same kind of smoothing as above, but for the distal end fo the fiber
    if np.any(phiShapeEmpirical[np.argmin(phiShapeEmpirical):]>0):
                
        last = np.where(phiShapeEmpirical[np.argmin(phiShapeEmpirical):]>0)[0][0]+np.argmin(phiShapeEmpirical)
        
        phiShapeEmpirical[last:] = 0
        
    else:
        last = np.where(np.abs(np.diff(phiShapeEmpirical))>1e-4)[0][-1]
        lasta = np.where(phiShapeEmpirical[last]+1e-4*np.arange(len(phiShapeEmpirical)-last)>0)[0][0]
        lasta += last

        phiShapeEmpirical[last:lasta] = 1e-4* np.arange(lasta-last)+ phiShapeEmpirical[last]
        phiShapeEmpirical[lasta:] = 0
    

    return xvals, phiShapeEmpirical

def FitPhiShape(fascIdx,distance,femDirectory):
    
    ''' 
    This function creates an interpolation object for the recording exposure
    '''

    phi = pd.read_excel(femDirectory+str(fascIdx)+'_BetterConductivity.xlsx')
    
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