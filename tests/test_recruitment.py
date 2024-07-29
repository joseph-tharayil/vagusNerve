import numpy as np
from vagusNerve.recruitment import *
import quantities as pq

def test_loadTitrationFactors():

    stimulusDirectory = {
    "myelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_HighConductivity.xlsx',
    "unmyelinated":'/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_HighConductivity.xlsx'
}

    m, u = loadTitrationFactors(stimulusDirectory)

    assert len(m)==50*39
    assert len(u)==50*39

def test_removeDuplicates():

    values = np.array([0,1,2,3,4,4,5,6])
    newvals = removeDuplicates(values)

    np.testing.assert_equal(newvals,np.array([0,1,2,3,4,5,6]))

def test_getCdf():

    titrationFactors = np.arange(10)
    fascIdx = 0

    midpts, cdf = getCdf(titrationFactors, fascIdx,removeJumps=False)

    np.testing.assert_equal(midpts,titrationFactors)
    np.testing.assert_equal(cdf,np.arange(1,11)/10)

    midpts, cdf = getCdf(titrationFactors, fascIdx,removeJumps=True)

    np.testing.assert_equal(midpts,titrationFactors)
    np.testing.assert_equal(cdf,np.arange(1,11)/10)

def test_jumpRemover():

    midpts = np.array([1,1.1,20,21])

    midptsX = jumpRemover(midpts,0)

    np.testing.assert_equal(midptsX,np.array([1,1.1]))
