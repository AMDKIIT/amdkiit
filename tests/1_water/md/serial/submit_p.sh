#!/bin/bash -l
#SBATCH --job-name=amdkiit_h2o
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#SBATCH --output=log.o%j
#SBATCH --partition=hm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --time=12:00:00
#SBATCH --export=NONE
#SBATCH --mail-type=NONE
#SBATCH --no-requeue

module load ohpc
export AMD=/home/msccp23/amdkiit-main/source/build/amdkiit.x

INPUT=input.yaml
OUTPUT=amdkiit.out
MYNP=1
mpirun -n ${MYNP} $AMD $INPUT > $OUTPUT
