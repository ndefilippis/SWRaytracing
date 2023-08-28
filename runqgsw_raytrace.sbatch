#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8GB
#SBATCH --time=8:00:00
#SBATCH --job-name=qgsw_raytracing
#SBATCH --output=slurm_%A_%a.out
#SBATCH --array=0-9
  
module purge
module load matlab/2020b

Ugs=(0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0)

Ug=${Ugs[$SLURM_ARRAY_TASK_ID]}

rundir=$SCRATCH/swqg_raytracing/run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID/.*}
mkdir -p $rundir
cd qg_flow_ray_trace
cp fulspec.m g2k.m k2g.m interpolate.m read_field.m write_field.m $rundir
cp grid_U.m interpolate_U.m $rundir
cp qgsw_raytrace.m $rundir
cd $rundir
mkdir data

matlab -nodisplay -r "qgsw_raytrace(256, 50, 2, 2000, 1200, $Ug, 3, 1); exit()" > run.log

exit

