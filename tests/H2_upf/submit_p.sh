#!/bin/bash -l
#SBATCH --job-name=amdkiit_h2o
#SBATCH -e ./%x.%j.err
#SBATCH -o ./%x.%j_log.o
#SBATCH --partition=hm
#SBATCH --nodes=1
#SBATCH --tasks-per-node=48
#SBATCH --time=00:15:00
#SBATCH --export=NONE
#SBATCH --mail-type=NONE
#SBATCH --no-requeue



module unload openmpi3/3.1.4
module unload gnu8/8.3.0

module load spack/0.17.1
source /home/apps/spack/share/spack/setup-env.sh
spack load intel-oneapi-compilers@2021.4.0
spack load intel-oneapi-mkl@2023.2.0

export AMD=/home/paramitag/amkdiit_intel/amdkiit-main/Github_upload_8jan24/source/build/amdkiit.x

INPUT=input.yaml
OUTPUT=${SLURM_JOB_NAME}.out
MYNP=48
mpirun -n ${MYNP} $AMD $INPUT > $OUTPUT
