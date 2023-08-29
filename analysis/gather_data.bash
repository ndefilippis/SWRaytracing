#!/bin/bash

slurm_job=run-$1
for i in {0..24}
do
  local_dir=job-$1/run-$i
  remote_dir=/scratch/nad9961/swqg_raytracing/run-$1_$i
  if [[ ! -e $local_dir ]]; then
    mkdir -p $local_dir
  fi
  scp nad9961@greene:$remote_dir/data/packet_*.bin $local_dir
  scp nad9961@greene:$remote_dir/run.log $local_dir
done
