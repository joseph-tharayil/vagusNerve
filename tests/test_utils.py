import numpy as np
from vagusNerve.utils import *
import quantities as pq

def test_velocity():

    d = np.array([20e-6])*pq.m
    v = getVelocities(d)
    np.testing.assert_almost_equal(v[0].item(),86.95)
    assert v[0].units == pq.m/pq.s

def test_scaling():

    d = (2/np.pi**.5)*pq.m
    fiberType = 1

    s = Scaling(d,fiberType)

    np.testing.assert_almost_equal(s.item(),1)

    assert s.units == pq.m/pq.ohm