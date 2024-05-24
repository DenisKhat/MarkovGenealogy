#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --job-name=Example
#SBATCH --mem=1024MB
#SBATCH --qos=privileged
#SBATCH --account=def-hpcg1709

echo $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR
cd Example

#  Put your job commands after this line
Rscript --no-save --no-restore slurm_methodB.R $SLURM_ARRAY_TASK_ID