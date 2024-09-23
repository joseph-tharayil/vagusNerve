# Apache-2.0

import numpy as np
import pandas as pd


import sys

from vagusNerve.runSim import runSim

def main(outputfolder, distanceIdx):

    stimulus = {'current':np.array([500])/173,
                'stimulusDirectory':{
                    "myelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationBetterConductivity_Standoff_HighConductivity_NoTime.xlsx',
                    "unmyelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/Titration_BetterSundt_HighConductivity.xlsx'
                }
               }

    recording = {'recordingCurrent':509e-6,
                 'recordingDirectory':'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/PhiConductivity_Bipolar_Corrected/'
            }

    runSim(outputfolder, distanceIdx, stimulus, recording)


if __name__=="__main__":
    
    
    outputfolder = sys.argv[1]
    distanceIdx = int(sys.argv[2])
    
    main(outputfolder,distanceIdx)
