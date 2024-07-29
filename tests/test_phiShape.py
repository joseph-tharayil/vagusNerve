import numpy as np
from vagusNerve.phiShape import *
import quantities as pq

# def test_phiShape():

#     recording = {'recordingCurrent':509e-6,
#                  'recordingDirectory':'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/PhiConductivity_Bipolar_Corrected/'
#             }

#     femDirectory = recording['recordingDirectory']

#     fascIdx = 0
#     distance = 0.06 * pq.m

#     fit = FitPhiShape(fascIdx,distance,femDirectory)

#     t = np.arange(100)*pq.s
#     velocity = [1*pq.m/pq.s]

#     p = PhiShape(velocity,t,fit)

#     assert p.units == pq.V
