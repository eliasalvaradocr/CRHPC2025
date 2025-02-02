#!/bin/bash
#SBATCH --job-name=collectives
#SBATCH --partition=nu
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:20:00
#SBATCH --nodes=1

module load gcc/11.1.0
module load mpich/4.1.1

CXX=mpic++

$CXX -o findpi findpi_solution_1.cpp

mpiexec -np 4 ./findpi 100000

rm findpi

$CXX -o findpi findpi_solution_2.cpp

mpiexec -np 4 ./findpi 100000

rm findpi

$CXX -o findpi findpi_solution_3.cpp

mpiexec -np 4 ./findpi 100000

rm findpi