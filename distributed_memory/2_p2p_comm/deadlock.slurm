#!/bin/bash
#SBATCH --job-name=helloworld_serial
#SBATCH --partition=nu
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:20:00
#SBATCH --nodes=1

module load gcc/11.1.0
module load mpich/4.1.1

CXX=mpic++

$CXX -o deadlock deadlock.cpp

mpirun -np 2 ./deadlock

# Get the PID of the last background process (MPI process)
#MPI_PID=$!

# Sleep for 5 seconds
#sleep 5

#echo "5 seconds passed"

# Kill the MPI process if it's still running
#kill $MPI_PID


#rm deadlock
