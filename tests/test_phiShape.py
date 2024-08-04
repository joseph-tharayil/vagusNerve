import numpy as np
from vagusNerve.phiShape import *
import quantities as pq

def test_removeZeroCrossings():

    testVector = np.array([-1,0,1,10,1,0,-10,-1,0,1])

    output = removeZeroCrossings(testVector)

    np.testing.assert_equal(output,np.array([0,0,1,10,1,0,-10,-1,0,0]))

def test_linearizeLeftSide():

    testVector = np.ones(7)
    cutoffPoint = 2
    slope = 1

    output = linearizeLeftSide(testVector, cutoffPoint, slope)

    expectedOutput = np.array([0,0,1,1,1,1,1])

    np.testing.assert_equal(output,expectedOutput)

def test_linearizeRightSide():

    testVector = -np.ones(7)
    cutoffPoint = 2
    slope = 1

    output = linearizeRightSide(testVector, cutoffPoint, slope)

    expectedOutput = np.array([-1,-1,-1,0,0,0,0])

    np.testing.assert_equal(output,expectedOutput)

def test_phiShape():

    recording = {'recordingCurrent':509e-6,
                 'recordingDirectory':'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/PhiConductivity_Bipolar_Corrected/'
            }

    femDirectory = recording['recordingDirectory']

    fascIdx = 0
    distance = 0.06 * pq.m

    fit = FitPhiShape(fascIdx,distance,femDirectory)

    t = np.arange(100)*pq.s
    velocity = [1*pq.m/pq.s]

    p = PhiShape(velocity,t,fit)

    np.testing.assert_equal(p.shape,(1,len(t)))
