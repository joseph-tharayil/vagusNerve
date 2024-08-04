from vagusNerve.runSim import *
import numpy as np
import quantities as pq

def test_getDiameterScalingOfCurrent():

    diameters = np.array([2e-6,2e-6])*pq.m
    time = np.array([0,1])*pq.s
    velocities = np.array([1,1])*pq.m/pq.s

    m, u = getDiameterScalingOfCurrent(diameters,time,velocities)

    assert m.units == pq.s**2/pq.ohm
    assert u.units == pq.s**2/pq.ohm

    assert len(m) == 2

def test_getDistance():

    recording = {}
    d = getDistance(0,recording)
    assert d.units == pq.m
    np.testing.assert_almost_equal(d.item(),0.06)

def test_getTime():

    t = getTime()
    assert t.units == pq.s

def test_getDiameters():

    d = getDiameters()
    assert d.units == pq.m

def test_getVariance():

    v = getVariance({})
    np.testing.assert_array_equal(v,np.array([0]))
