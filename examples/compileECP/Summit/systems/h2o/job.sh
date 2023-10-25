#!/bin/bash

#BSUB -P CSC304
#BSUB -W 01:10
#BSUB -nnodes 1
#BSUB -alloc_flags NVME
#BSUB -J uo2
#BSUB -o uo2o.%J
#BSUB -e uo2e.%J

module purge
module load cmake
module load cuda
module load gcc
module load netlib-lapack
module load openblas
module load magma
module load spectrum-mpi

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MEMBERWORK/csc304/Summit/metis-5.1.0/lib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MEMBERWORK/csc304/Summit/metis-5.1.0/build/Linux-x86_64/libmetis/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MEMBERWORK/csc304/Summit/magma/lib

# disable all the IBM optimized barriers and drop back to HCOLL or OMPI's barrier implementations
#export OMPI_MCA_coll_ibm_skip_barrier=true

# ROMIO has a hint for GPFS named IBM_largeblock_io which optimizes I/O with operations on large blocks
#export IBM_largeblock_io=true

# OpenMP: 1 thread per MPI rank
export OMP_NUM_THREADS=7

# run lmp
cd $MEMBERWORK/csc304/Summit/uo2/
jsrun -n1 -r1 -a1 -EOMP_NUM_THREADS=7 -g1 -c7 -brs $MEMBERWORK/csc304/Summit/lammps/src/lmp_serial < in.md
