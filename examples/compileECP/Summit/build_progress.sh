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

export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}

MY_PATH=$(pwd)

BML_LIB="${MY_PATH}/bml/install"
METIS_LIB="${MY_PATH}/metis-5.1.0/"

export MAGMA_ROOT=${OLCF_MAGMA_ROOT:="${OLCF_MAGMA_ROOT}"}
export PKG_CONFIG_PATH="$BML_LIB/lib/pkgconfig:$BML_LIB/lib64/pkgconfig"
export CMAKE_PREFIX_PATH=${CMAKE_PREFIX_PATH:="$BML_LIB:$BML_LIB/lib/:$BML_LIB/lib64/:$METIS_LIB/lib"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_MPI=${PROGRESS_MPI:=no}
export MAGMA_PATH=${MAGMA_PATH:="$OLCF_MAGMA_ROOT/"}
export CMAKE_INCLUDE_PATH=${CMAKE_INCLUDE_PATH:="$METIS_LIB/include"}
export CMAKE_LIBRARY_PATH=${CMAKE_LIBRARY_PATH:="$METIS_LIB/lib"}
export BLAS_VENDOR=${BLAS_VENDOR:=Generic}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/qmd-progress/install"}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export EXTRA_FCFLAGS=${EXTRA_FCFLAGS:=" -Wl,--copy-dt-needed-entries -I${MAGMA_PATH}/include -I${BML_LIB}/include/  -ffree-line-length-none -fopenmp -lpthread"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:="-Wl,--copy-dt-needed-entries -fopenmp  -L${BML_LIB}/lib64/ -lbml_fortran -lbml -L${MAGMA_PATH}/lib/ -lmagma  -L$CUDA_TOOLKIT_ROOT_DIR/lib64/ -lcublas -lcusolver -L${OLCF_CUDA_ROOT}/lib64/ -lcublas -lcudart  -fopenmp -fopenmp -lpthread  "}

cd qmd-progress; ./build.sh configure; cd build; make install; cd $MY_PATH
