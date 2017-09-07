cmake_minimum_required(VERSION 3.1)
project(latte Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/Modules)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CXX_FLAGS)

set(PRECISION DOUBLE CACHE STRING "The float precision")
string(TOUPPER ${PRECISION} PRECISION)

if((NOT ${PRECISION} STREQUAL SINGLE) AND (NOT ${PRECISION} STREQUAL DOUBLE))
  message(FATAL_ERROR "PRECISION needs to be either SINGLE or DOUBLE")
endif()

set(OPENMP FALSE CACHE BOOL "Use OpenMP")

if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

set(DO_MPI FALSE CACHE BOOL "Use MPI")
if(DO_MPI)
  find_package(MPI REQUIRED)
  if (NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Can not find suitable MPI library")
  else()
    message(STATUS "Found MPI library")
  endif()
  #set(MPI_Fortran_COMPILE_FLAGS "-DDO_MPI -I$ENV{MPI_ROOT}/include -I$ENV{MPI_ROOT}/lib")
  message(STATUS "MPI compile flags are ${MPI_Fortran_COMPILE_FLAGS}")
  message(STATUS "MPI includes are ${MPI_Fortran_INCLUDE_PATH}")
  message(STATUS "MPI libs are ${MPI_Fortran_LIBRARIES}")
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "Linking BLAS via ${BLAS_LIBRARIES}")
message(STATUS "Linking LAPACK via ${LAPACK_LIBRARIES}")
include(CheckFortranSourceCompiles)
set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS})
set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
check_Fortran_source_compiles("      program blas\n      call dgemm()\n      end program blas" BLAS_WORKS)
check_Fortran_source_compiles("      program lapack\n      call dsyev()\n      end program lapack" LAPACK_WORKS)
if(NOT BLAS_WORKS OR NOT LAPACK_WORKS)
  message(FATAL_ERROR "Can not compile against BLAS/LAPACK. Please see error log.")
endif()

set(XSYEV TRUE CACHE BOOL "Use DSYEV.")
set(XSYEVD FALSE CACHE BOOL "Use DSYEVD.")

if((NOT XSYEV AND NOT XSYEVD) OR (XSYEV AND XSYEVD))
  message(FATAL_ERROR "Please specify either XSYEV or XSYEVD")
endif()

option(DBCSR_OPT "Whether to use DBCSR" OFF)

option(GPUOPT "Whether to use the GPU" OFF)
if(GPUOPT)
  find_package(CUDA REQUIRED)
endif()

set(BML_BASE_DIR "/usr/local" CACHE STRING "The base path to the bml library")
find_library(BML_LIBS bml PATHS ${BML_BASE_DIR}/lib)
if(BML_LIBS)
  find_file(BML_HEADER_FOUND bml.mod PATHS ${BML_BASE_DIR}/include)
endif()
if(BML_LIBS AND BML_HEADER_FOUND)
  set(BML_FOUND TRUE)
  get_filename_component(BML_INCLUDE ${BML_HEADER_FOUND} DIRECTORY)
  message(STATUS "Found bml library: ${BML_LIBS}")
else()
  message(STATUS "Did not find bml library")
endif()

set(SOLVERS_DIR "/usr" CACHE STRING "The path to the qmmd-solvers library")
find_library(SOLVERS_LIBS qmmd-solvers PATHS ${SOLVERS_DIR}/lib)
if(SOLVERS_LIBS)
  find_file(SOLVERS_HEADER_FOUND qmmd-solvers.mod PATHS ${SOLVERS_DIR}/include)
endif()
if(SOLVERS_LIBS AND SOLVERS_HEADER_FOUND)
  set(SOLVERS_FOUND TRUE)
  get_filename_component(SOLVERS_INCLUDE ${SOLVERS_HEADER_FOUND} DIRECTORY)
  message(STATUS "Found qmmd-solvers library: ${SOLVERS_LIBS}")
else()
  message(STATUS "Did not find qmmd-solvers library")
endif()

string(TOLOWER latte_${PRECISION} EXECUTABLE_NAME)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMKAKE_SOURCE_DIR}/../src/Doxyfile.in Doxyfile)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
endif()

if(BML_FOUND)
  message(STATUS "Compiling bml dependent functions")
  set(BML_SOURCES
    ${CMKAKE_SOURCE_DIR}/../src/genX_bml_mod.F90
    ${CMKAKE_SOURCE_DIR}/../src/sp2bml_mod.F90
    ${CMKAKE_SOURCE_DIR}/../src/nonortho_bml_mod.F90
)
endif()

if(DBCSR_OPT STREQUAL ON)
  set(DBCSR_SOURCES 
    ${CMKAKE_SOURCE_DIR}/../src/dbcsr_var_mod.F90 
    ${CMKAKE_SOURCE_DIR}/../src/init_dbcsr.F90 
    ${CMKAKE_SOURCE_DIR}/../src/shutdown_dbcsr.F90
  )
endif()

add_executable(latte
  ${BML_SOURCES}
  ${DBCSR_SOURCES}
  ${CMKAKE_SOURCE_DIR}/../src/genX_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/allocatepure.F90
  ${CMKAKE_SOURCE_DIR}/../src/atomcharge.F90
  ${CMKAKE_SOURCE_DIR}/../src/conjgradient.F90
  ${CMKAKE_SOURCE_DIR}/../src/constants_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/controls_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/coulomb_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/deallocateall.F90
  ${CMKAKE_SOURCE_DIR}/../src/deallocatepure.F90
  ${CMKAKE_SOURCE_DIR}/../src/delta.F90
  ${CMKAKE_SOURCE_DIR}/../src/diag_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/diagdata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/energy_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/entropy.F90
  ${CMKAKE_SOURCE_DIR}/../src/fermi_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/fermidata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/fitting_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/force_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/gershgorin.F90
  ${CMKAKE_SOURCE_DIR}/../src/getdensity.F90
  ${CMKAKE_SOURCE_DIR}/../src/getmaxf.F90
  ${CMKAKE_SOURCE_DIR}/../src/getpressure.F90
  ${CMKAKE_SOURCE_DIR}/../src/getrho.F90
  ${CMKAKE_SOURCE_DIR}/../src/hbuild_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/homolumo_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/ifrestart.F90
  ${CMKAKE_SOURCE_DIR}/../src/initialize_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/io_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/kspace_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/kspacedata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/latte.F90
  ${CMKAKE_SOURCE_DIR}/../src/matrixio_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/md_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/mddata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/msrelax.F90
  ${CMKAKE_SOURCE_DIR}/../src/myprecision.F90
  ${CMKAKE_SOURCE_DIR}/../src/neblist_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/neblistdata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/nnz.F90
  ${CMKAKE_SOURCE_DIR}/../src/nonortho_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/nonorthodata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/norms.F90
  ${CMKAKE_SOURCE_DIR}/../src/openfiles_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/pairpot_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/pairpotdata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/panic_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/parallel_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/perturbation_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/plot_ppot.F90
  ${CMKAKE_SOURCE_DIR}/../src/plot_univ.F90
  ${CMKAKE_SOURCE_DIR}/../src/printsparse.F90
  ${CMKAKE_SOURCE_DIR}/../src/pulay_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/purearray.F90
  ${CMKAKE_SOURCE_DIR}/../src/random_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/relaxcommon.F90
  ${CMKAKE_SOURCE_DIR}/../src/resetprodhd.F90
  ${CMKAKE_SOURCE_DIR}/../src/restartarray.F90
  ${CMKAKE_SOURCE_DIR}/../src/scaling_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/scf_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/setuparray.F90
  ${CMKAKE_SOURCE_DIR}/../src/solvematcg.F90
  ${CMKAKE_SOURCE_DIR}/../src/solvematlapack.F90
  ${CMKAKE_SOURCE_DIR}/../src/sp2_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/sparsedata_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/sparsemath_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/sparsesp2_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/sparsity.F90
  ${CMKAKE_SOURCE_DIR}/../src/spinarray.F90
  ${CMKAKE_SOURCE_DIR}/../src/stdescent.F90
  ${CMKAKE_SOURCE_DIR}/../src/subgraph_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/subgraphsp2_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/timer_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/twonorm.F90
  ${CMKAKE_SOURCE_DIR}/../src/virialarray.F90
  ${CMKAKE_SOURCE_DIR}/../src/xbo_mod.F90
  ${CMKAKE_SOURCE_DIR}/../src/xbodata_mod.F90)

if(BML_HEADER_FOUND)
  include_directories(${BML_INCLUDE})
endif()

if(BML_LIBS)
  target_link_libraries(latte ${BML_LIBS})
endif()

if(SOLVERS_HEADER_FOUND)
  include_directories(${SOLVERS_INCLUDE})
endif()

if(SOLVERS_LIBS)
  target_link_libraries(latte ${SOLVERS_LIBS})
endif()

if(MPI_Fortran_INCLUDE_PATH)
  target_link_libraries(latte ${MPI_Fortran_INCLUDE_PATH})
endif()

if(MPI_Fortran_LIBRARIES)
  target_link_libraries(latte ${MPI_Fortran_LIBRARIES})
endif()

target_compile_definitions(latte
  PRIVATE
  ${PRECISION}PREC
  GPU${GPUOPT}
  DBCSR_${DBCSR_OPT}
  $<$<BOOL:${XSYEV}>:XSYEV>
  $<$<BOOL:${XSYEVD}>:XSYEVD>
  $<$<BOOL:${BML_FOUND}>:HAVE_BML>
  $<$<BOOL:${SOLVERS_FOUND}>:HAVE_SOLVERS>)

target_link_libraries(latte
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
  ${CUDA_cublas_LIBRARY}
  ${CUDA_cusparse_LIBRARY}
  ${CUDA_CUDART_LIBRARY})

set_target_properties(latte
  PROPERTIES
  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  OUTPUT_NAME ${EXECUTABLE_NAME})

if(OPENMP_FOUND)
  set_target_properties(latte
    PROPERTIES
    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
endif()

if(MPI_Fortran_FOUND)
  set_target_properties(latte
    PROPERTIES
    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
endif()

install(TARGETS latte
  RUNTIME DESTINATION bin)
