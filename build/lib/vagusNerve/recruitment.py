import numpy as np
import pandas as pd

from scipy.interpolate import interp1d

import quantities as pq

import sys


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


def loadTitrationFactors(stimulusDirectory):

        #### Loads and sorts titration factors from S4L. Sorting is justg to make sure that fibers and fascicles are in numerical order (ie, fiber 0-fiber50, fascicle0-fascicle39)

    titrationFactorsMeff = sortTitrationSpace(pd.read_excel(stimulusDirectory['myelinated'],index_col=0)).iloc[-1].values
    
    titrationFactorsUaff = sortTitrationSpace(pd.read_excel(stimulusDirectory['unmyelinated'],index_col=0)).iloc[-1].values

    return titrationFactorsMeff, titrationFactorsUaff

def removeDuplicates(midptsX):

    ### Removes points which have same titration factor. This removes vertical segments in the cdf, which can cause problems with interplation
    
    dupIdx =  np.where(np.diff(midptsX)==0)
    
    midptsX = np.delete(midptsX,dupIdx)

    return midptsX

def jumpRemover(midptsX, fascIdx):

    ### Finds and removes plateaus in the CDF. These occur because of a handful of fibers which have much higher thresholds than the others, which is not realistic

    diff = np.diff(midptsX/midptsX[0])
            
    jumpIdx = np.where(diff > 1.25)[0]
            
    if fascIdx != 35 and len(jumpIdx)>0:
        if len(jumpIdx)>1:
            jumpIdx = jumpIdx[0]
            
        end = len(midptsX)
        jumpRange = np.arange(jumpIdx+1,end)

        midptsX = np.delete(midptsX,jumpRange)

    return midptsX

def getCdf(titrationFac, fascIdx,removeJumps=True):

    ### Defines CDF of the titration curves

    midptsX = np.sort(titrationFac)
        
    midptsX = removeDuplicates(midptsX)

    if removeJumps:
    
        midptsX = jumpRemover(midptsX, fascIdx)

    cdfX = np.arange(1,len(midptsX)+1)/len(midptsX)

    midptsX = np.insert(midptsX,0,0)
    cdfX = np.insert(cdfX,0,0)

    return midptsX, cdfX


def interpolateTitrationFactors(titrationFac, current, diameters, d0, fascIdx,removeJumps=True):

    midptsX, cdfX = getCdf(titrationFac, fascIdx,removeJumps)
    
    interp = interp1d(midptsX,cdfX,bounds_error=False,fill_value=(0,1))

    thresholds = []

    for d in diameters:

        thresholds.append(interp(current*(d/d0)))

    return thresholds


def Recruitment(current,diameters, fascIdx,stimulusDirectory):
    
    d0Myelinated = 4e-6*pq.m
    d0Unmyelinated = 0.8e-6*pq.m

    titrationFactorsMeff, titrationFactorsUaff = loadTitrationFactors(stimulusDirectory)


    titrationFacM = np.array(titrationFactorsMeff[fascIdx*50:(fascIdx+1)*50]) # Selects fibers in fascicle
    titrationFacU = np.array(titrationFactorsUaff[fascIdx*50:(fascIdx+1)*50]) # Selects fibers in fascicle
        
    myelinated = interpolateTitrationFactors(titrationFacM, current, diameters, d0Myelinated, fascIdx,removeJumps=True)
    unmyelinated = interpolateTitrationFactors(titrationFacU, current, diameters, d0Unmyelinated, fascIdx,removeJumps=False)
    
    return [myelinated,unmyelinated]