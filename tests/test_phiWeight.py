# Apache-2.0

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
import numpy as np

def test_getPhiWeight():

    '''Tests that the function produces an output of the correct size'''

    current = np.arange(5)
    d = getDiameters()
    fascIdx = 0
    fascTypes = getFascicleTypes()

    distributionParams = {'maff':{'diameterParams':None, 'fiberTypeFractions':None},'meff':{'diameterParams':None, 'fiberTypeFractions':None}}


    stimulusDirectory = {
    "myelinated":'../../data/TitrationGoodConductivity_Standoff_Sideways_HighConductivity.xlsx',
    "unmyelinated":'../../data/TitrationGoodConductivity_Standoff_Sideways_Unmyelinated_HighConductivity.xlsx'
}

    phiWeight0, phiWeight1, phiWeight2, phiWeight3 = getPhiWeight(d, current,fascIdx,fascTypes, stimulusDirectory, distributionParams)

    np.testing.assert_array_equal(phiWeight0.shape,(len(current),len(d)))
    np.testing.assert_array_equal(phiWeight1.shape,(len(current),len(d)))
    np.testing.assert_array_equal(phiWeight2.shape,(len(current),len(d)))
    np.testing.assert_array_equal(phiWeight3.shape,(len(current),len(d)))