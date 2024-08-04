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

from vagusNerve.runSim import runSim

def main(outputfolder, distanceIdx):

    stimulus = {'current':np.arange(100,510,100)/28.6,
                'stimulusDirectory':{
                    "myelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_HighConductivity.xlsx',
                    "unmyelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_HighConductivity.xlsx'
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
