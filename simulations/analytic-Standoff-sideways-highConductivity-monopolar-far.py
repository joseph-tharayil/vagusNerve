# Apache-2.0

import numpy as np
import pandas as pd


import sys

from vagusNerve.runSim import runSim

def main(outputfolder, distanceIdx):

    stimulus = {'current':np.array([500])/28.6,
                'stimulusDirectory':{
                    "myelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_HighConductivity_NoTime.xlsx',
                    "unmyelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_HighConductivity.xlsx'
                }
               }

    recording = {'recordingCurrent':452e-6,
                 'recordingDirectory':'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/PhiConductivity_Monopolar_Far_Corrected/',
                 'distances': np.array([0.06,.05,0.01]) + 3 * 3e-3,
                 'cutoff':5e-6
            }

    runSim(outputfolder, distanceIdx, stimulus, recording)


if __name__=="__main__":
    
    
    outputfolder = sys.argv[1]
    distanceIdx = int(sys.argv[2])
    
    main(outputfolder,distanceIdx)
