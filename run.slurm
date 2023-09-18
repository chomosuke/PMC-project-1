#!/bin/bash

#SBATCH --job-name=ass1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=72
#SBATCH --time=00:00:5
#SBATCH --partition=cascade

echo 'Ohayo gazaimasu 72 threads'

module load Clang/13.0.1

make clean
make solution

OMP_NUM_THREADS=72 ./solution < grid10x10expensive.txt

make clean
