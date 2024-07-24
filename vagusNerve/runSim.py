import numpy as np
import pandas as pd

from mpi4py import MPI

import quantities as pq

from vagusNerve.phiWeight import *
from vagusNerve.utils import *
from vagusNerve.nerveSetup import *
from vagusNerve.phiShape import *

def runSim(outputfolder, distanceIdx, stimulus, recording):
    
   
    ## Selects fascicle and random seed for each available cpu
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    iteration = int(rank/39)
    fascIdx = int(rank % 39)
    #####

    current = stimulus['current']
    
    distances = [0.06,.05,0.01] # Stimulus-recording distance, in m
         
    
    distance = distances[distanceIdx]
    
    ### Loads action potential shapes in time
    ap = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShape20.xlsx') # Rat
    ap2 = pd.read_excel('/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Data/APShapePoint8.xlsx') # Sundt
    ####

    aps = [ap,ap2]

    nx=500000

    tmin=-3 # In s
    tmax=3 # In s
    tphi=np.arange(tmin,tmax,(tmax-tmin)/(nx-1))*pq.s


    recordingCurrent = recording['recordingCurrent'] # Current in the S4L recording simulation

    d0List = [20*1e-6,0.8*1e-6] # Diameters of myelinated and unmyelinated fibers used to calculate velocities

    velocities = [86.95,0.416] # Velcities for the above diamters
    
    
    fascTypes = getFascicleTypes(iteration) # Defines whether fasicle is on left or right side of nerve

    d = getDiameters(iteration)   

    stimulusDirectory = stimulus['stimulusDirectory']


    phiWeight0, phiWeight1, phiWeight2, phiWeight3 = getPhiWeight(d,current,fascIdx, fascTypes, stimulusDirectory) # Scaling factor for each diameter


    phi = [0,0,0,0]


    velocityList = getVelocities(d0List,velocities,d) # Gets velocity for each diameter


    names = ['myelinated','unmyelinated']

    recordingDirectory = recording['recordingDirectory']

    phiFunc = FitPhiShape(fascIdx,distance, recordingDirectory)# Defines an interpolation function for the recording exposure for the fasicle
    
    ### For each diameter, defines a shifted and scaled exposure function
    phiShape0 = PhiShape(velocityList[0],tphi,phiFunc)
    
    phiShape1 = PhiShape(velocityList[1],tphi,phiFunc)
    
#     np.save(outputfolder+'/'+str(iteration)+'/phis/'+str(distanceIdx)+'/rawShapes_'+str(fascIdx)+'.npy',[phiShape0,phiShape1])
   #############


### Scales exposure functions

    scaling0 = (Scaling(d,0)* (tphi[1]-tphi[0])/velocityList[0])
    scaling1 = (Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])
    
    np.save(outputfolder+'/'+str(iteration)+'/diameters/scaling_'+str(fascIdx)+'.npy',[scaling0,scaling1])
    
    scaling00 = phiWeight0.T*scaling0[:,np.newaxis]
    scaling10 = phiWeight1.T*scaling0[:,np.newaxis]
    scaling01 = phiWeight2.T*scaling1[:,np.newaxis]
    scaling11 = phiWeight3.T*scaling1[:,np.newaxis]
    
    np.save(outputfolder+'/'+str(iteration)+'/maff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling00)
    
    np.save(outputfolder+'/'+str(iteration)+'/meff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling10)
    
    np.save(outputfolder+'/'+str(iteration)+'/uaff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling01)
    
    np.save(outputfolder+'/'+str(iteration)+'/ueff/'+str(distanceIdx)+'/scaling_'+str(fascIdx)+'.npy',scaling11)
    
    phi[0] += np.matmul(phiShape0.T,scaling00)
    
    phi[1] += np.matmul(phiShape0.T,scaling10)
    
    phi[2] += np.matmul(phiShape1.T,phiWeight2.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])
    
    phi[3] += np.matmul(phiShape1.T,phiWeight3.T*(Scaling(d,1)* (tphi[1]-tphi[0])/velocityList[1])[:,np.newaxis])


    phi = np.array(phi)
    
    np.save(outputfolder+'/'+str(iteration)+'/'+'phis/'+str(distanceIdx)+'/'+str(fascIdx)+'.npy',phi)
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
