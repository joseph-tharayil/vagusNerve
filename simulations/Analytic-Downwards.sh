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
#SBATCH --output=EEG_0_CoordsV.out
#SBATCH --error=EEG_0_CoordsV.err
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --array=0-1

module purge
module load archive/2022-02 py-mpi4py

source ~/probevenv/bin/activate

filename='/gpfs/bbp.cscs.ch/project/proj85/scratch/vagusNerve/Analytic_allActive_squareTerm_downward'

rm -r $filename
mkdir $filename


mkdir $filename/0
mkdir $filename/0/diameters
mkdir $filename/0/fascicles
mkdir $filename/0/recruitment
mkdir $filename/0/phis

mkdir $filename/0/maff
mkdir $filename/0/maff/0
mkdir $filename/0/maff/1

cp -r $filename/0/maff $filename/0/meff
cp -r $filename/0/maff $filename/0/uaff
cp -r $filename/0/maff $filename/0/ueff

srun -n 39 python analytic-Standoff-highConductivity.py $filename $SLURM_ARRAY_TASK_ID
wait
srun -n 4 python combineOldMethod.py $filename $SLURM_ARRAY_TASK_ID

