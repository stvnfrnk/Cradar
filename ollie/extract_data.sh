#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --array=1-24
#SBATCH --output=/home/ollie/sfranke/extract_%A_%a.out
#SBATCH -p smp

##Enlarge the stacksize, just to be on the safe side.
ulimit -s unlimited

echo "SLURM_JOBID: "$SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: "$SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: "$SLURM_ARRAY_JOB_ID

chan=$(ls /work/ollie/sfranke/Data/$1/UWB/ | grep 'chan' | awk -v line=$SLURM_ARRAY_TASK_ID '{if (NR==line) print $0}')
srun tar xvfj /work/ollie/sfranke/Data/$1/UWB/${chan} -C /work/ollie/sfranke/Data/

