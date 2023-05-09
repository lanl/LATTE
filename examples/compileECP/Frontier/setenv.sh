module purge  
module load PrgEnv-gnu
module load magma/2.6.1
module load cmake
module load openblas
module load rocm/5.2.0
export LD_LIBRARY_PATH=${ROCM_PATH}/lib:${LD_LIBRARY_PATH}

