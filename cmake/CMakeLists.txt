cmake_minimum_required(VERSION 3.1)
project(latte Fortran)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)
  #release comes with -O3 by default
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE)
endif(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_Fortran_FLAGS)

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

option(DBCSR_OPT "Whether to use DBCSR" OFF)

option(GPUOPT "Whether to use the GPU" OFF)
if(GPUOPT)
  find_package(CUDA REQUIRED)
endif()

string(TOLOWER latte_${PRECISION} EXECUTABLE_NAME)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_SOURCE_DIR}/../src/Doxyfile.in Doxyfile)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
endif()

set(MYMODULES myprecision.f90 constants_mod.f90 setuparray.f90 coulombarray.f90
diagarray.f90 fermicommon.f90 mdarray.f90 neblistarray.f90
nonoarray.f90 ppotarray.f90 purearray.f90 relaxcommon.f90 restartarray.f90
sparsearray.f90 spinarray.f90 univarray.f90 virialarray.f90 xboarray.f90
kspacearray.f90 timer_mod.f90 subgraph.f90 matrixio.f90 sparsemath.f90
sparsesp2.f90 subgraphsp2.f90 homolumo.f90 openfiles_mod.f90 kernelparser_mod.f90
latteparser_latte_mod.f90 mixer_mod.f90 genXprogress.f90 latte_lib.f90)

set(DBCSRMODS dbcsr_var_mod.f90 init_dbcsr.f90 shutdown_dbcsr.f90)

set(MYSUBS allocatediag.f90 allocatenono.f90 allocatepure.f90 bodirect.f90
deallocatediag.f90 deallocatepure.f90 deorthomyrho.f90 deallocatenono.f90
diagmyh.f90 entropy.f90 fermiexpans.f90 genX.f90  latte.f90 orthomyH.f90
 orthomyrho.f90  pulay.f90 pulay_sp.f90 solvematcg.f90 solvematlapack.f90
sp2T.f90
sp2fermi.f90 sp2fermi_init.f90 sp2pure.f90 sp2pure_sparse.f90 kdiagmyh.f90
 kbodirect.f90 sp2gap.f90 sp2gap_setup.f90 readrestart.f90 initrng.f90
parawrite.f90 parafileopen.f90 wrtcfgs.f90 wrtrestart.f90 pbc.f90
addqdep.f90 allocatecoulomb.f90 allocatenebarrays.f90  allocatexbo.f90 
atomcharge.f90 avepress.f90 avetemp.f90 bldnewH.f90 bldspinH.f90
conjgradient.f90 coulomb_ewald.f90 coulomb_oldskool.f90 coulomb_rspace.f90
coultailcoef.f90 deallocateall.f90 deallocatecoulomb.f90
deallocatenebarrays.f90  deallocatexbo.f90 fcoulnono.f90 fcoulnono_sp.f90
fermiallocate.f90 fermideallocate.f90
fspinnono.f90 fspinnono_sp.f90  gershgorin.f90 getbndfil.f90
getcoule.f90 getdensity.f90 getdeltaq.f90 getdeltaspin.f90 gethdim.f90
getke.f90 getmaxf.f90 getmdf.f90 getpressure.f90 getrho.f90 getspinE.f90 
 ifrestart.f90 initcoulomb.f90 initialv.f90 initshockcomp.f90 msrelax.f90
neblists.f90 noelec.f90 norms.f90 nptrescale.f90  nvtrescale.f90
pairpot.f90 panic.f90 plot_ppot.f90 plot_univ.f90 printsparse.f90
 propchempot_xbo.f90 propspins_xbo.f90 qconsistency.f90 readcontrols.f90
readcr.f90 readmdcontroller.f90 readppot.f90  readtb.f90 resetprodhd.f90
rhozero.f90 shockcomp.f90  spinrhodirect.f90 stdescent.f90
summary.f90 tbmd.f90 toteng.f90 univscaling_function.f90 univtailcoef.f90 
vdwtailcoef.f90 velverlet.f90  xbo.f90 am.f90 bm.f90 factorial.f90
 slmmp.f90 tlmmp.f90 wignerd.f90 dfda.f90 dfdb.f90 dfdr.f90
dunivscaling_function.f90 dslmmpda.f90 dslmmpdb.f90 dtlmmpda.f90
dtlmmpdb.f90 dwignerddb.f90 gradH.f90 kbldnewh.f90  kgetrho.f90 kgradH.f90 
bldnewHS_sp.f90 gradH_sp.f90 univscaling.f90 dunivscaling.f90 kgetdos.f90
qneutral.f90 shiftH.f90 getmatindlist.f90 nvtlangevin.f90 nvtandersen.f90
nvtNH.f90 sp2pure_sparse_parallel.f90 sp2pure_sparse_parallel_simple.f90
sp2pure_subgraph_parallel.f90 gaussrn.f90 nnz.f90 allocatesubgraph.f90
deallocatesubgraph.f90 dosfit.f90 mofit.f90 getrespf.f90 gendiag.f90
mofit_plato.f90 ppfit.f90 pairpot_noneb.f90 allfit.f90
readppottab.f90 pairpottab.f90 hugrescale.f90 avesforhug.f90 getdipole.f90
fittingoutput.f90 tabtest.f90 kgenX.f90 korthomyH.f90 kdeorthomyrho.f90
kfcoulnono.f90 kpulay.f90 getforce.f90 masses2symbols.f90 setuptbmd.f90
wrtrestartlib.f90 readrestartlib.f90) 

foreach(_m MYMODULES MYSUBS DBCSR)
  set(${_m}_SOURCES)
  foreach(_file ${${_m}})
    list(APPEND ${_m}_SOURCES ${CMAKE_SOURCE_DIR}/../src/${_file})
  endforeach()
endforeach()

if(NOT DBCSR_OPT)
  set(DBCSR_SOURCES)
endif()

add_library(latte ${MYMODULES_SOURCES} ${MYSUBS_SOURCES})
target_compile_options(latte PRIVATE -cpp)

  #if(SOLVERS_HEADER_FOUND)
  #  include_directories(${SOLVERS_INCLUDE})
  #endif()
  #
  #if(SOLVERS_LIBS)
  #  target_link_libraries(latte ${SOLVERS_LIBS})
  #endif()
  #
  #if(MPI_Fortran_INCLUDE_PATH)
  #  target_link_libraries(latte ${MPI_Fortran_INCLUDE_PATH})
  #endif()
  #
  #if(MPI_Fortran_LIBRARIES)
  #  target_link_libraries(latte ${MPI_Fortran_LIBRARIES})
  #endif()
  #
target_compile_definitions(latte
  PRIVATE
    XSYEV
    ${PRECISION}PREC
    GPU${GPUOPT}
    DBCSR_${DBCSR_OPT}
    MPI_${DO_MPI}
    $<$<BOOL:${DO_MPI}>:DO_MPI_ON>
)
  #
  #target_link_libraries(latte
  #  ${LAPACK_LIBRARIES}
  #  ${BLAS_LIBRARIES}
  #  ${CUDA_cublas_LIBRARY}
  #  ${CUDA_cusparse_LIBRARY}
  #  ${CUDA_CUDART_LIBRARY})
  #
  #set_target_properties(latte
  #  PROPERTIES
  #  RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}
  #  OUTPUT_NAME ${EXECUTABLE_NAME})
  #
  #if(OPENMP_FOUND)
  #  set_target_properties(latte
  #    PROPERTIES
  #    COMPILE_FLAGS ${OpenMP_Fortran_FLAGS}
  #    LINK_FLAGS ${OpenMP_Fortran_FLAGS})
  #endif()
  #
  #if(MPI_Fortran_FOUND)
  #  set_target_properties(latte
  #    PROPERTIES
  #    COMPILE_FLAGS ${MPI_Fortran_COMPILE_FLAGS})
  #endif()
  #
  #install(TARGETS latte
  #  RUNTIME DESTINATION bin)
