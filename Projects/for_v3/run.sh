#!/bin/bash -l

#SBATCH -A snic2019-3-432
#SBATCH -p core -n 8
#SBATCH -t 00:00:05
#SBATCH -J pcmc

#SBATCH -e error_file.e
#SBATCH -o output_file.o

mpirun -np 8 ./pcmc

