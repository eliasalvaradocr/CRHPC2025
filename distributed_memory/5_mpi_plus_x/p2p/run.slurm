#!/bin/bash
#SBATCH --job-name=collectives
#SBATCH --partition=nu-wide
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:20:00
#SBATCH --nodes=2

module --force purge
module load gcc/11.1.0
module load mpich/4.1.1
module load cmake/3.26.0
module load kokkos/openmp

mpiexec -np 2 ./build/p2p