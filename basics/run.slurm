#!/bin/bash

#SBATCH --job-name=<program>
#SBATCH --partition=kura
#SBATCH --output=%x-stdout_%j
#SBATCH --error=%x-stderr_%j
#SBATCH --time=0-00:20:00

#SBATCH --nodes=1


if [ $SLURM_JOB_PARTITION == 'kura' ]; then
  . /opt/Modules/3.2.10/init/sh
fi

module purge
# Load needed modules
if [[ $SLURM_JOB_PARTITION == 'nukwa' || $SLURM_JOB_PARTITION == 'nukwa-v100' || $SLURM_JOB_PARTITION == 'nukwa-debug' ]]; then
  module load cuda/11.1
  module load kokkos/cuda
  nvidia-smi -L
else
  module load kokkos/openmp
fi
module load gcc/11.1.0
module load cmake/3.26.0

