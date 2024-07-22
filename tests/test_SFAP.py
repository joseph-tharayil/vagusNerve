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

import quantities as pq

def test_SFAP():

    ap = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShape20.xlsx')

    ap2 = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapePoint8.xlsx')

    aps = [ap,ap2]

    nx=1000000

    tmin=-3
    tmax=3
    tphi=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))#*pq.s
    
    d0List = np.array([20*1e-6,0.8e-6,1*1e-6])#*pq.m
    
    velocities = np.array([86.95,0.413])#*pq.m/pq.s
    
    ds = np.array([20e-6,0.8e-6,0.4e-6])#*pq.m

    segmentLength= 1e-6#*pq.m

    Vs = getVs(aps,tphi)
        
    k = np.linspace(-1,1-(2.0/np.shape(Vs)[1]),np.shape(Vs)[1])

    for fiberType in range(1):
        
        d = 1e-6#*pq.m #ds[0]

        if fiberType == 0:
            deff = d #.32*d
        else:
            deff = d

        velocityList = getVelocities(d0List,velocities,d)
        
        distances = [.10757,0.0176]#*pq.m #distances[0]
        
        if fiberType == 0:
            distance = distances[0]
        else:
            distance = distances[1]

        v = velocityList[fiberType]

        
        V = Vs[fiberType]

        der = np.diff(V,n=2)/((tphi[1]-tphi[0])**2)

        phiFunc = FitPhiShape(0,distance, '/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/FasciclePhiBigStandard/')

        phiShape0 = PhiShape(velocityList,tphi,phiFunc)[0]

        cv = fftconvolve(der*Scaling(deff,d0List[fiberType],fiberType)* (tphi[1]-tphi[0])/v,phi,mode='same') 
        
