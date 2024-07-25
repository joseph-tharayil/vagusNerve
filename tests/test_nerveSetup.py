import numpy as np
from vagusNerve.nerveSetup import *

def test_getFasciclePositions():

    '''Tests that the function returns an array of shape number_fascicles x 3'''

    pos = getFasciclePositions()

    np.testing.assert_array_equal(pos.shape,(39,3))

def test_getFascicleTypes():

    '''Tests that the function produces a reasonable number of fascicles on each side'''

    fascTypes = getFascicleTypes()

    numberPerSide = np.sum(fascTypes)

    assert int(numberPerSide)==numberPerSide

    assert numberPerSide > 10 and numberPerSide < 30