#!/bin/bash

rm -r build
rm -r install

# Set BML Library location
MY_PATH=`pwd`
export BML_DIR=${MY_PATH}/bml/install
echo $BML_DIR
# Configuring PROGRESS with OpenMP
export MAGMA_PATH=${MAGMA_PATH:=${OLCF_MAGMA_ROOT}}
export CC=${CC:=gcc}
export FC=${FC:=gfortran}
export CXX=${CXX:=g++}
export PROGRESS_OPENMP=${PROGRESS_OPENMP:=yes}
export INSTALL_DIR=${INSTALL_DIR:="${MY_PATH}/qmd-progress/install"}
export PROGRESS_GRAPHLIB=${PROGRESS_GRAPHLIB:=no}
export PROGRESS_TESTING=${PROGRESS_TESTING:=yes}
export CMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE:=RelWithDebInfo}
export PROGRESS_EXAMPLES=${PROGRESS_EXAMPLES:=yes}
export PROGRESS_BENCHMARKS=${PROGRESS_BENCHMARKS:=yes}
export PKG_CONFIG_PATH=${BML_LIB}/lib64/pkgconfig
export EXTRA_FCFLAGS=""
export EXTRA_LINK_FLAGS=""
cd qmd-progress; ./build.sh configure; cd build; make install; cd $MY_PATH
