#!/bin/bash

# Make sure all the paths are correct

rm -r build
rm -r install

MY_PATH=$(pwd)

export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export BLAS_VENDOR=${BLAS_VENDOR:=OpenBLAS}
export BML_OPENMP=${BML_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/bml/install"}
export BML_MAGMA=${BML_MAGMA:=yes}
export MAGMA_ROOT=${OLCF_MAGMA_ROOT}
export BML_ROCSOLVER=${BML_ROCSOLVER:=yes}
export BML_TESTING=${BML_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=Release}
export CMAKE_Fortran_FLAGS=${CMAKE_Fortran_FLAGS:="-g -ffree-form -ffree-line-length-200"}
export CMAKE_C_FLAGS=${CMAKE_C_FLAGS:="-g -I${HIP_PATH}/include"}
export EXTRA_LINK_FLAGS=${EXTRA_LINK_FLAGS:=""}

if [[ -z "$CMAKE_PREFIX_PATH" ]]; then
    export CMAKE_PREFIX_PATH=${ROCM_PATH}
else
    export CMAKE_PREFIX_PATH="${ROCM_PATH};${CMAKE_PREFIX_PATH}"
fi

cd bml; ./build.sh configure; cd build; make -j; make install ; cd $MY_PATH

