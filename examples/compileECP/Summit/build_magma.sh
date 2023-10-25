!/bin/bash
# Compile MAGMA
MY_PATH=$(pwd)
export FPIC="-fPIC"
export CFLAGS='"-O3 ${FPIC} -DNDEBUG -DADD_  -DUPCASE -Wshadow"'
export FFLAGS='"-O3 ${FPIC}  -DADD_  -DUPCASE  -DNDEBUG"'
export F90FLAGS='"-O3 ${FPIC}  -DADD_  -DUPCASE    -DNDEBUG"'
export NVCCFLAGS='"-O3 -DNDEBUG -DADD_  -DUPCASE -Xcompiler -fPIC -std=c++11"'
export LDFLAGS='"-ffast-math -fPIC -fopenmp"'
export LIB='"-fopenmp -lpthread -lstdc++ -lm -lgfortran -lcublas -lcusparse -lcudart -lcudadevrt"'
export LIBDIR='"-L${CUDA_PATH}/lib64 "'
export INC='"-I${CUDA_INCLUDES}/ -I${MY_PATH}/magma/include "'
export OPTIONS=${OPTIONS:=" VERBOSE=1 CC=cc CXX=CC FORT=ftn CFLAGS=${CFLAGS} FFLAGS=${FFLAGS} F90FLAGS=${F90FLAGS} \
NVCCFLAGS=${NVCCFLAGS} LDFLAGS=${LDFLAGS} LIB=${LIB} LIBDIR=${LIBDIR} INC=${INC}"}

echo "make" $OPTIONS lib -j > ./magma/compile.sh
cd magma; cp ./make.inc-examples/make.inc.power9-essl make.inc
make clean; sh ./compile.sh

