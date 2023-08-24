#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4GB
#SBATCH --time=48:00:00
#SBATCH --job-name=sqgp1
#SBATCH --output=slurm_%j.out
  
module load matlab/2017a

rundir=$SCRATCH/sqgp1/run-${SLURM_JOB_ID/.*}
mkdir -p $rundir
cp sqgp1.m $rundir
cp runsqgp1.m $rundir
cd $rundir

matlab -nodisplay -r "runsqgp1; exit()" > run.log 


exit

