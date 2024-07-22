import numpy as np
from mpi4py import MPI
import sys

def main(outputfolder,distanceIdx):

    rank = MPI.COMM_WORLD.Get_rank()

    iteration = 0

    typeIdx = rank

    names = ['maff','meff','uaff','ueff']

    signals = [0,0,0,0]

    for fascIdx in range(39):

        signals[typeIdx] += np.load(outputfolder+'/'+str(iteration)+'/'+names[typeIdx]+'/'+str(distanceIdx)+'/signals_'+str(fascIdx)+'.npy')

    np.save(outputfolder+'/'+str(iteration)+'/'+names[typeIdx]+'/'+str(distanceIdx)+'/signals.npy',signals[typeIdx])


if __name__=='__main__':
    
    outputfolder = sys.argv[1]
    distanceIdx = int(sys.argv[2])
    
    main(outputfolder,distanceIdx)
    
    