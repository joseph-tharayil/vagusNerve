# Apache-2.0

import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
import sys

from vagusNerve.runSim import runSim

def runSim_wrapper(fascIdx, stim, rec,outputfolder):
    params = {'maff':{'diameterParams':None, 'fiberTypeFractions':None},'meff':{'diameterParams':None, 'fiberTypeFractions':None}}
    return runSim(fascIdx, stim, rec, params, 2000,outputfolder)  # Pass correct arguments

def main(outputfolder):

    stimulus = {'current':np.array([100,200,300,400,500])/173,
                'stimulusDirectory':{
                    "myelinated":'../../titrationTestVertical.xlsx'
                }
               }

    recording = {'recordingCurrent':509e-6,
                 'recordingDirectory':'../../SimResults/'
            }

    numcores = mp.cpu_count()
    with mp.Pool(numcores-18) as p:
        signals = p.starmap(runSim_wrapper, [(i, stimulus, recording,outputfolder) for i in np.arange(39)])

    np.save(outputfolder+'/results.npy',signals)

if __name__=="__main__":
    
    
    outputfolder = 'downwardsElectrode/'
    
    main(outputfolder)
