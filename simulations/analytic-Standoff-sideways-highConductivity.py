# Apache-2.0

import numpy as np
import pandas as pd


import sys

from vagusNerve.runSim import runSim

def main(outputfolder, distanceIdx):

    stimulus = {'current':np.arange(100,510,100)/27.7,
                'stimulusDirectory':{
                    "myelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_BetterHighConductivity_NoTime.xlsx',
                    "unmyelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_BetterHighConductivity.xlsx'
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
