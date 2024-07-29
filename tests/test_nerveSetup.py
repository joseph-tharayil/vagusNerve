import numpy as np
from vagusNerve.nerveSetup import *
from vagusNerve.utils import getDiameters

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

def test_prob():

    '''
    Tests that probability function produces output of the same size as the input vector
    '''

    d = getDiameters()
    maffProb = 1
    p = MaffProb(d,maffProb)

    assert len(p) == len(d)

def test_gamma_dist():

    '''
    Tetsts that gamma distribution function produces accurate results
    '''

    x = 4
    k = 9
    theta = 0.5

    np.testing.assert_almost_equal(gammaDist(x,k,theta), 0.28,decimal=3)

def test_fiberTypeFractions():

    '''
    Tests that we get reasonable fractions for each fiber type
    '''

    fascIdx = 0
    fascTypes = [1]

    maffFrac, meffFrac, ueffFrac, uaffFrac = getFiberTypeFractions(fascIdx, fascTypes)

    assert maffFrac < .2
    assert maffFrac > 0
    assert meffFrac > 0.05
    assert meffFrac < .2
    assert meffFrac > maffFrac
    assert ueffFrac > 0
    assert ueffFrac < .3

    assert uaffFrac == 1-(maffFrac+meffFrac+ueffFrac)

    fascTypes = [0]

    maffFrac, meffFrac, ueffFrac, uaffFrac = getFiberTypeFractions(fascIdx, fascTypes)

    assert maffFrac < .25
    assert maffFrac > .10
    assert meffFrac > 0.
    assert meffFrac < .05
    assert maffFrac > meffFrac
    assert ueffFrac > .05
    assert ueffFrac < .3

    assert uaffFrac == 1-(maffFrac+meffFrac+ueffFrac)

def test_fiberNumbers():

    fascTypes = getFascicleTypes()

    numbers = []
    for i in range(39):
        numbers.append(getFibersPerFascicle(i,fascTypes))

    assert np.sum(numbers)>250000
    assert np.sum(numbers)<375000
    

    