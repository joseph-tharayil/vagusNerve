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
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --array=0-1

# Apache-2.0

filename='/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/results/finalResults/Analytic_downward'

if [ $SLURM_ARRAY_TASK_ID -eq 0 ]

then

echo $SLURM_ARRAY_TASK_ID

rm -r $filename
mkdir $filename


mkdir $filename/diameters
mkdir $filename/fascicles
mkdir $filename/recruitment
mkdir $filename/phis

mkdir $filename/phis/0
mkdir $filename/phis/1
mkdir $filename/phis/2

mkdir $filename/maff
mkdir $filename/maff/0
mkdir $filename/maff/1
mkdir $filename/maff/2

mkdir $filename/meff
mkdir $filename/meff/0
mkdir $filename/meff/1
mkdir $filename/meff/2

mkdir $filename/uaff
mkdir $filename/uaff/0
mkdir $filename/uaff/1
mkdir $filename/uaff/2

mkdir $filename/ueff
mkdir $filename/ueff/0
mkdir $filename/ueff/1
mkdir $filename/ueff/2

fi

srun -n 39 python analytic-Standoff-downwards-highConductivity.py $filename $SLURM_ARRAY_TASK_ID
wait
srun -n 4 python combineOldMethod.py $filename $SLURM_ARRAY_TASK_ID

