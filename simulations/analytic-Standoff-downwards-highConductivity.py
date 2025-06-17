# Apache-2.0

import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
import sys

from vagusNerve.runSim import runSim

def runSim_wrapper(fascIdx, stim, rec):
    params = {'maff':{'diameterParams':None, 'fiberTypeFractions':None},'meff':{'diameterParams':None, 'fiberTypeFractions':None}}
    return runSim(0, stim, rec, fascIdx,params, 2000)  # Pass correct arguments

def main(outputfolder):

    stimulus = {'current':np.array([100,200,300,400,500])*10/173,
                'stimulusDirectory':{
                    "myelinated":r'D:/vagusNerve/FullPipeline/titrationTestVertical.xlsx'
                }
               }

    recording = {'recordingCurrent':509e-6,
                 'recordingDirectory':r'C:/Users/tharayil/Desktop/SimResults/'
            }

    numcores = mp.cpu_count()
    with mp.Pool(numcores-4) as p:
        signals = p.starmap(runSim_wrapper, [(i, stimulus, recording) for i in np.arange(39)])

    np.save(outputfolder+'/results.npy',signals)


if __name__=="__main__":
    
    
    outputfolder = 'downwardsElectrode/'
    
    main(outputfolder)
