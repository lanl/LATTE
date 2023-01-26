cmake_minimum_required(VERSION 3.9)
project(latteqeq VERSION 1.2.2 LANGUAGES Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/Modules)

include(GNUInstallDirs)

set(PRECISION DOUBLE CACHE STRING "The float precision")
set_property(CACHE PRECISION PROPERTY STRINGS DOUBLE SINGLE)
string(TOUPPER ${PRECISION} PRECISION)

if((NOT ${PRECISION} STREQUAL SINGLE) AND (NOT ${PRECISION} STREQUAL DOUBLE))
  message(FATAL_ERROR "PRECISION needs to be either SINGLE or DOUBLE")
endif()

option(OPENMP "Use OpenMP" OFF)
if(OPENMP)
  find_package(OpenMP REQUIRED)
endif()

option(DO_MPI "Use MPI" OFF)
if(DO_MPI)
  find_package(MPI REQUIRED)
  if (NOT MPI_Fortran_FOUND)
    message(FATAL_ERROR "Can not find suitable MPI library")
  endif()
endif()

find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

message(STATUS "Linking BLAS via ${BLAS_LIBRARIES}")
message(STATUS "Linking LAPACK via ${LAPACK_LIBRARIES}")
#include(CheckFortranSourceCompiles)
#set(CMAKE_REQUIRED_FLAGS ${OpenMP_Fortran_FLAGS})
#set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${BLAS_LIBRARIES})
#check_Fortran_source_compiles("      program blas\n      call dgemm()\n      end program blas" BLAS_WORKS)
#check_Fortran_source_compiles("      program lapack\n      call dsyev()\n      end program lapack" LAPACK_WORKS)
#if(NOT BLAS_WORKS OR NOT LAPACK_WORKS)
#  message(FATAL_ERROR "Can not compile against BLAS/LAPACK. Please see error log.")
#endif()

set(MYMODULES latteqeq_lib.F90  qeq_parameters.F90 openfiles_mod.F90 kernelparser_mod.F90)

set(MYSUBS  Ewald_k_Space.F90 CoulombMatrix.F90
	#LoadBondIntegralParameters_H.F90 BondIntegral.F90
 KK0Res.F90  Ewald_Real_Space.F90
 AtomicDensityMatrix.F90
 Invert.F90
 #LoadPairPotParameters.F90  Eig.F90
 #SlaterKosterForce.F90  Get_dH.F90  Get_dS_tmp.F90
 #Canon_DM_PRT.F90  Get_dH_save.F90  Kernel_Fermi_NoDiag.F90
 #LoadBondIntegralParameters_S.F90  Build_S.F90  Build_H0.F90
 #ScaleTail.F90 PairPotential_Forces.F90
 nearestneighborlist.F90  Initiate.F90  masses2symbols.F90
 latteqeq_lib.F90
 openfiles_mod.F90
 kernelparser_mod.F90 
 latteqeq_parser.F90
 latteqeq_c_bind.F90  qeq_parameters.F90)

foreach(_m MYMODULES MYSUBS)
  set(${_m}_SOURCES)
  foreach(_file ${${_m}})
    list(APPEND ${_m}_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../src/${_file})
  endforeach()
endforeach()

add_library(latteqeq ${MYMODULES_SOURCES} ${MYSUBS_SOURCES}
	${CMAKE_CURRENT_SOURCE_DIR}/../src/latteqeq_lib.F90
	${CMAKE_CURRENT_SOURCE_DIR}/../src/latteqeq_c_bind.F90
)
target_compile_options(latteqeq PRIVATE -cpp)
install(TARGETS latteqeq EXPORT LATTEQEQ_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT LATTEQEQ_Targets FILE LATTEQEQ_Targets.cmake NAMESPACE LATTEQEQ:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTEQEQ)
include(CMakePackageConfigHelpers)
configure_file(LATTEQEQConfig.cmakein ${CMAKE_CURRENT_BINARY_DIR}/LATTEQEQConfig.cmake @ONLY)
write_basic_package_version_file("LATTEQEQConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/LATTEQEQConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/LATTEQEQConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTEQEQ)
add_library(LATTEQEQ::latteqeq ALIAS latteqeq)
target_link_libraries(latteqeq PRIVATE
  ${PROGRESS_LIBRARIES}
  ${BML_LIBRARIES}
  ${CUDA_cublas_LIBRARY}
  ${CUDA_cusparse_LIBRARY}
  ${CUDA_CUDART_LIBRARY}
  ${MPI_Fortran_LIBRARIES}
  PUBLIC
  ${LAPACK_LIBRARIES}
  ${BLAS_LIBRARIES}
)
if(OPENMP)
  target_link_libraries(latteqeq PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_include_directories(latteqeq PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${PROGRESS_INCLUDE_DIRS} ${BML_INCLUDE_DIRS})

#only "ON" and "OFF" support
foreach(_def GPUOPT DO_MPI PROGRESS)
  if(${${_def}})
    set(${_def} ON)
  else()
    set(${_def} OFF)
  endif()
endforeach()

target_compile_definitions(latteqeq PRIVATE XSYEV ${PRECISION}PREC GPU${GPUOPT} MPI_${DO_MPI} PROGRESS${PROGRESS})

if(DO_MPI)
  set(LATTEQEQ LATTEQEQ_MPI_${PRECISION})
else()
  set(LATTEQEQ LATTEQEQ_${PRECISION})
endif()

#add_executable(${LATTEQEQ} ${CMAKE_CURRENT_SOURCE_DIR}/../src/DFTB_POLFF.F90)
#target_link_libraries(${LATTEQEQ} latteqeq ${PROGRESS_LIBRARIES})
#target_include_directories(${LATTEQEQ} PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${PROGRESS_INCLUDE_DIRS} ${BML_INCLUDE_DIRS})
#target_compile_definitions(${LATTEQEQ} PRIVATE XSYEV ${PRECISION}PREC GPU${GPUOPT} MPI_${DO_MPI} PROGRESS${PROGRESS})
#install(TARGETS ${LATTEQEQ} DESTINATION ${CMAKE_INSTALL_BINDIR})
