#Forces the build system to have cmake from a certain version above!
cmake_minimum_required(VERSION 3.9)

project(latte VERSION 1.2.2 LANGUAGES Fortran)

#Enables the -DPROGRESS=<on/off> option?
option(PROGRESS "Use progress library" OFF)

#Since progress depends on BML and bml is writtend in C, then:
if(PROGRESS)
  enable_language(C CXX)
endif()

if(PROGRESS)
  include(FindPkgConfig)
  #Looks for the installation of progress in CMAKE_MODULE_PATH or CMAKE_PREFIX_PATH or whatever else cmake needs!
  find_package(PROGRESS CONFIG QUIET)
  #Checks if the PROGRESS modules are indeed available
  pkg_check_modules(PROGRESS REQUIRED progress)
  message(STATUS "Found progress: ${PROGRESS_LDFLAGS}")
  #Appending the necesary flags into LINK_LIBRARIES
  list(APPEND LINK_LIBRARIES ${PROGRESS_LDFLAGS})
  list(APPEND PROGRESS_LIBRARIES ${PROGRESS_LDFLAGS})

  #Same thing for BML!
  find_package(BML CONFIG QUIET)
  pkg_check_modules(BML REQUIRED bml)
  list(APPEND LINK_LIBRARIES ${BML_LDFLAGS})
  list(APPEND BML_LIBRARIES ${BML_LDFLAGS})
  message(STATUS "Found bml: ${BML_LDFLAGS}")
else()
  find_package(BLAS REQUIRED)
  find_package(LAPACK REQUIRED)
  message(STATUS "Linking BLAS via ${BLAS_LIBRARIES}")
  message(STATUS "Linking LAPACK via ${LAPACK_LIBRARIES}")
endif()

list(APPEND LINK_LIBRARIES "-fopenmp -lm -Wl,--copy-dt-needed-entries")

#If build type (release,debub,etc) is not specified, then set it as Release
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif()

#Cmake modules/macros are in a subdirectory to keep this file cleaner
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

#In case installation path is not specify the following will set /usr/bin, /usr/local/bin, ect
include(GNUInstallDirs)

set(PRECISION DOUBLE CACHE STRING "The float precision")
set_property(CACHE PRECISION PROPERTY STRINGS DOUBLE SINGLE)
string(TOUPPER ${PRECISION} PRECISION)

#Raise and error if presicion is not set
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


