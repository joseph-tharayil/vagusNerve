import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

import quantities as pq

import sys

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *


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



def Recruitment(current,diameters, fascIdx,stimulusDirectory):
    
    d0Myelinated = 4e-6*pq.m
    d0Unmyelinated = 0.8e-6*pq.m
    
   
    #### Loads and sorts titration factors from S4L. Sorting is justg to make sure that fibers and fascicles are in numerical order (ie, fiber 0-fiber50, fascicle0-fascicle39)

    titrationFactorsMeff = sortTitrationSpace(pd.read_excel(stimulusDirectory['myelinated'],index_col=0)).iloc[-1].values
    
    titrationFactorsUaff = sortTitrationSpace(pd.read_excel(stimulusDirectory['unmyelinated'],index_col=0)).iloc[-1].values

    ####


    titrationFactors = [titrationFactorsMeff, titrationFactorsUaff]


    for j in [0,1]: # Myelinated and unmyelinated, respectively

        titrationFac = np.array(titrationFactors[j][fascIdx*50:(fascIdx+1)*50]) # Selects fibers in fascicle
        
        
        midptsX = np.sort(titrationFac)
        
        dupIdx =  np.where(np.diff(midptsX)==0)
        
        midptsX = np.delete(midptsX,dupIdx)
        
        cdfX = np.arange(0,len(midptsX))/len(midptsX)
        
        diff = np.diff(midptsX/midptsX[0])
                
        jumpIdx = np.where(diff > 1.25)[0]
                
        if fascIdx != 35 and j == 0 and len(jumpIdx)>0:
            if len(jumpIdx)>1:
                jumpIdx = jumpIdx[0]
                
            end = len(midptsX)
            jumpRange = np.arange(jumpIdx,end)

            midpts2 = np.delete(midptsX,jumpRange)
            cdf2 = np.arange(0,len(midpts2))/len(midpts2)
            
#            cdfX = cdf2
#            midptsX = midpts2
            
        
        if j == 0:
            
            midpts = midptsX
            cdf = cdfX

        if j == 1:
            
            midptsU = midptsX
            cdfU = cdfX

### Defines CDF of the titration curves
    
    interp = interp1d(midpts,cdf,bounds_error=False,fill_value=(0,1))

    interpU = interp1d(midptsU,cdfU,bounds_error=False,fill_value=(0,1))
##############
        
    myelinated = []
    unmyelinated = []


    for j, d in enumerate(diameters):

        myelinated.append(interp(current*(d/d0Myelinated)))

        unmyelinated.append(interpU(current*d/d0Unmyelinated))

    
    return [myelinated,unmyelinated]