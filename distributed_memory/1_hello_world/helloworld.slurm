#!/bin/bash
#SBATCH --job-name=helloworld_serial
#SBATCH --partition=nu
#SBATCH --output=stdout-%x_%j
#SBATCH --error=stderr-%x_%j
#SBATCH --time=0-00:05:00
#SBATCH --nodes=1

module load gcc/11.1.0

CXX=/opt/compilers/gcc-11.1.0/bin/g++

$CXX -o helloworld helloworld.cpp

./helloworld

rm helloworld
