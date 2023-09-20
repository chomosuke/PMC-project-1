#!/bin/bash

#SBATCH --job-name=ass1
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:00:30
#SBATCH --partition=cascade

module load Clang/13.0.1

make clean
make solution

for OMP_NUM_THREADS in 1 2
do
	echo $OMP_NUM_THREADS threads
	for INPUT in $(find . -type f -name "grid*")
	do
		echo $INPUT
		./solution < $INPUT
	done
done
