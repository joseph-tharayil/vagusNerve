# Apache-2.0

import numpy as np
import pandas as pd


from scipy.stats import norm
from scipy.io import loadmat

from scipy.optimize import leastsq
from scipy.optimize import least_squares
from scipy.io import loadmat
from scipy.interpolate import interp1d
from scipy.stats import norm
import multiprocessing as mp
from scipy.fft import fft, ifft, fftshift,ifftshift
from scipy.signal import fftconvolve, sosfilt, butter

from vagusNerve.nerveSetup import *
from vagusNerve.utils import *
from vagusNerve.phiShape import *
from vagusNerve.runSim import *

import quantities as pq
from scipy.signal import find_peaks     

def test_SFAP():

    '''
    Calculates single-fiber action potential for a myelinated fiber with diameter 7.5 um. Tests that all of the units are correct, and that the features of the SFAP match those of the Sim4Life simulation
    '''

    aps = loadActionPotentialShapes()

    tphi=getTime()
    
    segmentLength= 1e-6*pq.m

    Vs = getVs(aps,tphi)
        
    fiberType = 0
    
    d = 7.5e-6*pq.m 

    deff = d 

    velocityList = getVelocities(d)
    
    distances = [.10757,0.0176]*pq.m
    
    distance = distances[0]

    v = velocityList[fiberType]

    V = Vs[fiberType]

    der = np.diff(V,n=2)/((tphi[1]-tphi[0])**2)

    assert der.units == pq.V/pq.s**2

    density_length = der*Scaling(deff,fiberType)*(1/v**2)

    assert density_length.units == pq.A/pq.m

    density_area = der*Scaling(deff,fiberType)*(1/v**2)/(np.pi*deff)

    assert density_area.units == pq.A/pq.m**2

    current = der*Scaling(d,fiberType)*(1/v**2)*segmentLength

    assert current.units == pq.A

    phiFunc = FitPhiShape(0,distance,'../../data/FasciclePhiBigStandard/')

    phiShape0 = PhiShape(velocityList,tphi,phiFunc)

    phi = phiShape0[fiberType]

    cv = fftconvolve(der*Scaling(deff,fiberType)* (tphi[1]-tphi[0])/v,phi,mode='same') 

    
    time = tphi[1:-1]*1e3

    cap = -cv/(726e-6)

    np.save('cap.npy',cap)
    
    positive_peaks = find_peaks(cap,height=0.25e-6)

    ppeakTimes = time[positive_peaks[0]]

    assert len(ppeakTimes)==2
    assert ppeakTimes[0] > 3 and ppeakTimes[0] < 4
    assert ppeakTimes[1] > 3 and ppeakTimes[1] < 4
    assert positive_peaks[1]['peak_heights'][0] > 2e-6 and positive_peaks[1]['peak_heights'][0] < 2.5e-6
    assert positive_peaks[1]['peak_heights'][1] > .2e-6 and positive_peaks[1]['peak_heights'][1] < .5e-6

    negative_peaks = find_peaks(-cap,height=0.25e-6)

    npeakTimes = time[negative_peaks[0]]

    assert len(npeakTimes)==2
    assert npeakTimes[0] < ppeakTimes[0]
    assert npeakTimes[1] > ppeakTimes[0] and npeakTimes[1] < ppeakTimes[1]
    assert negative_peaks[1]['peak_heights'][0] > .4e-6 and negative_peaks[1]['peak_heights'][0] < .55e-6
    assert negative_peaks[1]['peak_heights'][1] > 1.5e-6 and negative_peaks[1]['peak_heights'][1] < 2e-6
    
