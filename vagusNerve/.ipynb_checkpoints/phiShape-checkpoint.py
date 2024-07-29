import numpy as np
import pandas as pd


from scipy.interpolate import interp1d
from scipy.signal import fftconvolve, butter, sosfilt

import quantities as pq

def removeZeroCrossings(phiShapeEmpirical):

    '''
    If the phi curve crosses the zero line on the far ends of the nerve, sets the potential lateral to the zero crossing to zero
    '''

    if np.any(phiShapeEmpirical[:np.argmax(phiShapeEmpirical)]<0): # If the potential is negative at the end of the fiber, sets the potential to 0
        
        first = np.where(phiShapeEmpirical[:np.argmax(phiShapeEmpirical)]<0)[0][-1] 

        phiShapeEmpirical[:first+1] = 0

    
    if np.any(phiShapeEmpirical[np.argmin(phiShapeEmpirical):]>0): # Does the same for the other end of the fiber
                
        last = np.where(phiShapeEmpirical[np.argmin(phiShapeEmpirical):]>0)[0][0]+np.argmin(phiShapeEmpirical)
        
        phiShapeEmpirical[last:] = 0

    return phiShapeEmpirical

def linearizeRightSide(phiShapeEmpirical, cutoffPoint, slope):

    '''
    Given a point cutoffPoint on the phi curve, linearizes the region from cutoffPoint to the zeroCrossing, and sets points more lateral to the zero crossing to 0
    '''

    xIdx = np.arange(0,len(phiShapeEmpirical)-cutoffPoint)

    linearizedValues = phiShapeEmpirical[cutoffPoint]+ slope * xIdx

    phiShapeEmpirical[cutoffPoint:] = linearizedValues

    zeroCrossing = np.argmin(np.abs(phiShapeEmpirical[cutoffPoint:])) + cutoffPoint

    phiShapeEmpirical[zeroCrossing:] = 0 

    return phiShapeEmpirical

def linearizeLeftSide(phiShapeEmpirical, cutoffPoint, slope):

    '''
    Given a point cutoffPoint on the phi curve, linearizes the region from cutoffPoint to the zeroCrossing, and sets points more lateral to the zero crossing to 0
    '''

    xIdx = np.arange(-cutoffPoint,0)

    linearizedValues = phiShapeEmpirical[cutoffPoint]+ slope * xIdx

    phiShapeEmpirical[:cutoffPoint] = linearizedValues

    zeroCrossing = np.argmin(np.abs(phiShapeEmpirical[:cutoffPoint]))

    phiShapeEmpirical[:zeroCrossing] = 0 

    return phiShapeEmpirical

def smoothPhiShape(phiShapeEmpirical):

    phiShapeEmpirical = removeZeroCrossings(phiShapeEmpirical)

    ####### Makes sure that the potential at the proximal end of the fiber goes all the way to zero

    slopeToLinearize = 1e-4

    derivative = np.diff(phiShapeEmpirical)

    if phiShapeEmpirical[0] != 0: # If the potential does not go all the way to 0 by the end of the fiber, forces it to zero
        
        first = np.where(derivative>slopeToLinearize)[0][0] # Based on derivative of function, selects point after which not to change values

        phiShapeEmpirical = linearizeLeftSide(phiShapeEmpirical, first, slopeToLinearize)
        
    ############
   
    #### Does the same kind of smoothing as above, but for the distal end of the fiber
    if phiShapeEmpirical[-1] != 0:
        
        last = np.where(derivative>slopeToLinearize)[0][-1]
        
        phiShapeEmpirical = linearizeRightSide(phiShapeEmpirical, last, slopeToLinearize)
        
    return phiShapeEmpirical
   
def editPhiShape(phi,distance):
    
    ''' 
    This function takes the recording exposure curve from S4L, shifts it to match the desired distance from stimulus to recording, and smooths it
    '''

    xPos = phi.iloc[:,0]*pq.m
    
    xvals = xPos+distance -xPos[np.argmax(phi.iloc[:,1].values)] # Shift to match desired distance

    phiShapeEmpirical = (phi.iloc[:,1].values-np.mean(phi.iloc[:,1]))

    phiShapeEmpirical = smoothPhiShape(phiShapeEmpirical)

   ######## 
    

    return xvals, phiShapeEmpirical*pq.V

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