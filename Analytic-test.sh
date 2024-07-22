#!/bin/bash -l
#SBATCH --job-name="EEG_2_CoordsV"
#SBATCH --partition=prod
#SBATCH --nodes=39
##SBATCH -C clx
#SBATCH --cpus-per-task=2
#SBATCH --time=24:00:00
##SBATCH --mail-type=ALL
#SBATCH --account=proj85
#SBATCH --no-requeue
#SBATCH --output=EEG_1_CoordsV.out
#SBATCH --error=EEG_1_CoordsV.err
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --array=0-2

spack env activate bluerecording-dev
source ~/bluerecording-dev/bin/activate

filename='/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/results/finalResults/Analytic-test'

echo $SLURM_ARRAY_TASK_ID

if [ $SLURM_ARRAY_TASK_ID -eq 0 ]

then

echo $SLURM_ARRAY_TASK_ID

rm -r $filename
mkdir $filename


mkdir $filename/0
mkdir $filename/0/diameters
mkdir $filename/0/fascicles
mkdir $filename/0/recruitment
mkdir $filename/0/phis

mkdir $filename/0/phis/0
mkdir $filename/0/phis/1
mkdir $filename/0/phis/2

mkdir $filename/0/maff
mkdir $filename/0/maff/0
mkdir $filename/0/maff/1
mkdir $filename/0/maff/2

mkdir $filename/0/meff
mkdir $filename/0/meff/0
mkdir $filename/0/meff/1
mkdir $filename/0/meff/2

mkdir $filename/0/uaff
mkdir $filename/0/uaff/0
mkdir $filename/0/uaff/1
mkdir $filename/0/uaff/2

mkdir $filename/0/ueff
mkdir $filename/0/ueff/0
mkdir $filename/0/ueff/1
mkdir $filename/0/ueff/2

fi

wait

srun -n 39 python analytic-Standoff-sideways-highConductivity.py $filename $SLURM_ARRAY_TASK_ID
wait
srun -n 4 python combineOldMethod.py $filename $SLURM_ARRAY_TASK_ID

