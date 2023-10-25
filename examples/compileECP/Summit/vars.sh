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
