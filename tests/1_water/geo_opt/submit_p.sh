#!/bin/bash -l
#SBATCH --job-name=amdkiit_h2o
#SBATCH -o ./%x.%j.out
#SBATCH -e ./%x.%j.err
#SBATCH --output=log.o%j
#SBATCH --partition=hm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --time=00:15:00
#SBATCH --export=NONE
#SBATCH --mail-type=NONE
#SBATCH --no-requeue

module load openmpi/openmpi_4.1.2
module load compiler/gcc/8.3.0 

export AMD=/home/msccp23/amdkiit-main/source/build/amdkiit.x

INPUT=input.yaml
OUTPUT=amdkiit.out
MYNP=48
mpirun -n ${MYNP} $AMD $INPUT > $OUTPUT
