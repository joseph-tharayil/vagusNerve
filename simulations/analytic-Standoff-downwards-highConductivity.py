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
from scipy.signal import fftconvolve, butter, sosfilt

from scipy.stats import rv_histogram
from mpi4py import MPI

from math import gamma

from scipy.optimize import curve_fit

import quantities as pq

import sys

    
def main(outputfolder,distanceIdx):
    
   
    ## Selects fascicle and random seed for each available cpu
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    iteration = int(rank/39)
    fascIdx = int(rank % 39)
    #####
    
    distances = [0.06,.05] # Stimulus-recording distance, in m
    
    #### Random generators
    diameterGenerator = np.random.default_rng(seed=fascIdx)# Diameters depend on fascicle
    latencyGenerator = np.random.default_rng(seed=iteration) # Latency depends on trial index
    #####
       
    
    distance = distances[distanceIdx]
    
    ### Loads action potential shapes in time
    ap = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShape20.xlsx') # Rat
    ap2 = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapePoint8.xlsx') # Sundt
    ####

    aps = [ap,ap2]

    nx=500000

    tmin=-3 # In s
    tmax=3 # In s
    tphi=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))

    current = np.linspace(100,500,5)/173 # Currents from 100 to 1000 uA, scaled by the current in the S4L titration simulation
    

    recordingCurrent = 509e-6 # Current in the S4L recording simulation

    d0List = [20*1e-6,.8*1e-6] # Diameters of myelinated and unmyelinated fibers used to calculate velocities

    velocities = [86.95,0.416] # Velcities for the above diamters
    
    
    fascTypes = getFascicleTypes(iteration) # Defines whether fasicle is on left or right side of nerve

    d = getDiameters()   


    phiWeight0, phiWeight1, phiWeight2, phiWeight3 = getPhiWeight(d,current,fascIdx, fascTypes) # Scaling factor for each diameter


    phi = [0,0,0,0]


    velocityList = getVelocities(d0List,velocities,d) # Gets velocity for each diameter


    names = ['myelinated','unmyelinated']


    phiFunc = FitPhiShape(fascIdx,distance)# Defines an interpolation function for the recording exposure for the fasicle
    
    ### For each diameter, defines a shifted and scaled exposure function
    phiShape0 = PhiShape(velocityList[0],tphi,phiFunc)
    
    phiShape1 = PhiShape(velocityList[1],tphi,phiFunc)
   #############


### Scales exposure functions

    scaling0 = (Scaling(d,0)* (tphi[1]-tphi[0])/velocityList[0])
    
    np.save(outputfolder+'/'+str(iteration)+'/diameters/scaling_'+str(fascIdx)+'.npy',scaling0)
    
    scaling00 = phiWeight0.T*scaling0[:,np.newaxis]
    scaling10 = phiWeight1.T*scaling0[:,np.newaxis]
    
    np.save(outputfolder+'/'+str(iteration)+'/maff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling00)
    
    np.save(outputfolder+'/'+str(iteration)+'/meff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling10)
    
    phi[0] += np.matmul(phiShape0.T,scaling00)
    
    phi[1] += np.matmul(phiShape0.T,scaling10)
    
    phi[2] += np.matmul(phiShape1.T,phiWeight2.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])
    
    phi[3] += np.matmul(phiShape1.T,phiWeight3.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])


    phi = np.array(phi)
    
    np.save(outputfolder+'/'+str(0)+'/'+'phis/'+str(fascIdx)+'.npy',phi)
################


    signal = 0

    Vs = getVs(aps,tphi) # Interpolates action potential shape in time

    k = np.linspace(-1,1-(2.0/np.shape(Vs)[1]),np.shape(Vs)[1])

    signals = []

    for i, V in enumerate(Vs):  

        der = np.diff(V,n=2)/((tphi[1]-tphi[0])**2) # Second derivative of action potential shape

        cv = []
        
        for j in range(len(current)):
            
            c = fftconvolve(der,phi[i,:,j],mode='same') # Convolves second derivative with exposure
            
            cv.append(c)

        cv = np.array(cv)

        signals.append(cv)

    signals = np.array(signals)
    
    signals /= recordingCurrent


    names = ['maff','meff','uaff','ueff']

    for typeIdx in range(len(names)):

        np.save(outputfolder+'/'+str(iteration)+'/'+names[typeIdx]+'/'+str(distanceIdx)+'/signals_'+str(fascIdx)+'.npy',signals[typeIdx])


if __name__=="__main__":
    
    
    outputfolder = sys.argv[1]
    distanceIdx = int(sys.argv[2])
    
    main(outputfolder,distanceIdx)
