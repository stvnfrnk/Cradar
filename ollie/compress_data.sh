#!/bin/bash

#SBATCH --time=12:00:00
#SBATCH --output=/home/ollie/sfranke/compress_%A_%a.out
#SBATCH -p smp

##Enlarge the stacksize, just to be on the safe side.
ulimit -s unlimited

echo "SLURM_JOBID: "$SLURM_JOBID
echo "SLURM_ARRAY_TASK_ID: "$SLURM_ARRAY_TASK_ID
echo "SLURM_ARRAY_JOB_ID: "$SLURM_ARRAY_JOB_ID

DIR_IN="/work/ollie/sfranke/Data"
DIR_OUT="/work/ollie/sfranke/Data_compressed"

srun tar -zcvf ${DIR_OUT}/${1}-dms1/UWB/chan${2}.tar.gz ${DIR_IN}/${1}/chan${2}

