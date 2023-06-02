#!/bin/bash
#SBATCH -A CSC304
#SBATCH -J uo2
#SBATCH -o %x-%j.out
#SBATCH -t 00:30:00
#SBATCH -p batch
#SBATCH -N 1
export OMP_NUM_THREADS=8
module load PrgEnv-gnu
module load magma/2.6.1
module load cmake
module load openblas
module load rocm/5.2.0 

$MEMBERWORK/csc304/Frontier1/lammps/src/lmp_serial < in.md | tee out.log 
