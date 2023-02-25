cmake_minimum_required(VERSION 3.9)

option(PROGRESS "Use progress library" OFF)

project(latte VERSION 1.2.2 LANGUAGES Fortran)
if(PROGRESS)
  enable_language(C CXX)
endif()

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)

# Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

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

option(PROGRESS "Use progress library" OFF)
if(PROGRESS)
  find_package(PROGRESS CONFIG REQUIRED)
  if(TARGET PROGRESS::progress)
    set(PROGRESS_LIBRARIES PROGRESS::progress)
  endif()
  message(STATUS "Found PROGRESS: ${PROGRESS_LIBRARIES} ${BML_LIBRARIES}")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)

  message(STATUS "Linking BLAS via ${BLAS_LIBRARIES}")
  message(STATUS "Linking LAPACK via ${LAPACK_LIBRARIES}")
endif()  

option(DBCSR_OPT "Whether to use DBCSR" OFF)

option(GPUOPT "Whether to use the GPU" OFF)
if(GPUOPT)
  find_package(CUDA REQUIRED)
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/Doxyfile.in Doxyfile)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
endif()

if(DO_MPI)
  set(LATTE LATTE_MPI_${PRECISION})
else()
  set(LATTE LATTE_${PRECISION})
endif()

add_subdirectory(src)

#install library and executable
install(TARGETS ${LATTE} DESTINATION ${CMAKE_INSTALL_BINDIR})

install(TARGETS latte EXPORT LATTE_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
                                           ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT LATTE_Targets FILE LATTE_Targets.cmake NAMESPACE LATTE:: 
                             DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTE)

include(CMakePackageConfigHelpers)
configure_file(cmake/LATTEConfig.cmakein ${CMAKE_CURRENT_BINARY_DIR}/LATTEConfig.cmake @ONLY)
write_basic_package_version_file("LATTEConfigVersion.cmake" VERSION ${PROJECT_VERSION}
                                                            COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/LATTEConfig.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/LATTEConfigVersion.cmake"
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTE)

add_library(LATTE::latte ALIAS latte)
