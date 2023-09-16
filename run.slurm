#!/bin/bash

#SBATCH --job-name=ass1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:00:10
#SBATCH --partition=cascade

echo 'Ohayo gazaimasu'

module load OpenMPI/4.1.4
module load Clang/13.0.1

make clean
make solution

OMP_NUM_THREADS=4 ./solution < grid10x10expensive.txt

make clean
