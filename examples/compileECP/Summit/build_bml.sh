#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

module purge
module load cmake
module load cuda
module load gcc
module load netlib-lapack
module load openblas
module load magma
module load spectrum-mpi

MY_PATH=$(pwd)

export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}

export BLAS_VENDOR=${BLAS_VENDOR:=Generic}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/bml/install"}
export BML_TESTING=${BML_TESTING:=no}
export BML_COMPLEX=${BML_COMPLEX:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}

export MAGMA_PATH=${MAGMA_PATH:="$OLCF_MAGMA_ROOT/magma"}
export MAGMA_ROOT=${MAGMA_ROOT:="${OLCF_MAGMA_ROOT}"}
export BML_MAGMA=${BML_MAGMA:=yes}
export BML_CUSOLVER=${BML_CUSOLVER:=yes}

export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-fopenmp -lpthread  -L${OLCF_CUDA_ROOT}/lib64/ -lcublas -lcudart -ffixed-line-length-512"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-fopenmp -lpthread -L${OLCF_CUDA_ROOT}/lib64/ -lcublas -lcudart -ffixed-line-length-512"}



cd bml; ./build.sh configure; cd build; make -j; make install cd $MY_PATH


