#!/bin/bash
#SBATCH --job-name=random_walk
#SBATCH --partition=nukwa
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:20:00
#SBATCH --nodes=1

module --force purge
module load gcc/11.1.0
module load mpich/4.1.1
module load vtk/9.4.0
module load cuda/11.8
module load cmake/3.26.0
module load kokkos/cuda

./build/random_walk

tar -czf particles_random_walk.tar.gz *vtp*