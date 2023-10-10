#!/bin/bash
#PBS -N 32H2O
#PBS -q small
#PBS -l nodes=2:ppn=20
#PBS -j oe
cd $PBS_O_WORKDIR
export I_MPI_FABRICS=shm:dapl
export I_MPI_MPD_TMPDIR=/scratch/kritama
source /usr/mpi/gcc/openmpi-1.10.5a1/bin/mpivars.sh
AMDEXE=/home/kritama/AMDKIIT_05_09_2023/paramita/source/build/amdkiit.x   #path to the executable

MYNP=40
mpirun -machinefile $PBS_NODEFILE -n ${MYNP} $AMDEXE input.yaml > amdkiit.out

