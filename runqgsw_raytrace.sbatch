#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --job-name=qgsw_raytracing
#SBATCH --output=slurm_%A_%a.out
#SBATCH --array=1-20
  
module purge
module load matlab/2020b

config=parameters.txt

w0=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $2}' $config)
Ug=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $3}' $config)
f=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $4}' $config)
Cg=$(awk -v ArrayTaskID=$SLURM_ARRAY_TASK_ID '$1==ArrayTaskID {print $5}' $config)

rundir=$SCRATCH/swqg_raytracing/run-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID/.*}
mkdir -p $rundir
cd qg_flow_ray_trace
cp fulspec.m g2k.m k2g.m interpolate.m read_field.m write_field.m $rundir
cp grid_U.m interpolate_U.m $rundir
cp qgsw_raytrace.m $rundir
cd $rundir
mkdir data

matlab -nodisplay -r "qgsw_raytrace(256, 50, $w0, 6000, 1000, $Ug, $f, $Cg); exit()" > run.log

exit

