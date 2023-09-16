#!/bin/bash

#SBATCH --job-name=ass1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:00:30
#SBATCH --partition=cascade

echo 'Ohayo gazaimasu'

export OMP_NUM_THREADS=8

module load OpenMPI/4.1.4
module load Clang/13.0.1

make clean
make solution

./solution < grid10x10expensive.txt
