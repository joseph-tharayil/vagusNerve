# Apache-2.0

import numpy as np
import pandas as pd
import multiprocessing as mp
from functools import partial
import sys

from vagusNerve.runSim import runSim

def runSim_wrapper(fascIdx, stim, rec):
    params = {'maff':{'diameterParams':None, 'fiberTypeFractions':None},'meff':{'diameterParams':None, 'fiberTypeFractions':None}}
    return runSim(fascIdx, stim, rec, params, 2000)  # Pass correct arguments

def main(outputfolder):

    stimulus = {'current':np.array([500])/420.4,
                'stimulusDirectory': {
                    "myelinated": '../../titrationTestHorizontal.xlsx',
                }
                }

    recording = {'recordingCurrent':273e-6,
                 'recordingDirectory':'../../SimResults/PhiConductivity_Small_Otherside_20240215/',
                 'cutoff':1e-6
            }

    numcores = mp.cpu_count()
    with mp.Pool(12) as p:
        signals = p.starmap(runSim_wrapper, [(i, stimulus, recording) for i in np.arange(39)])

    np.save(outputfolder + '/results.npy', signals)


if __name__=="__main__":
    
    
    outputfolder = 'smallRec_otherSide'
    
    main(outputfolder)
