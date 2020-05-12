cmake_minimum_required(VERSION 3.9)
project(latte VERSION 1.2.2 LANGUAGES Fortran)

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

option(PROGRESS "Use progress library" OFF)
if(PROGRESS)
  if(NOT TARGET PROGRESS::progress)
    find_package(PROGRESS CONFIG QUIET)
  endif()
  if(TARGET PROGRESS::progress)
    set(PROGRESS_LIBRARIES PROGRESS::progress)
  else()
    find_package(BML MODULE REQUIRED)
    find_package(PROGRESS MODULE REQUIRED)
  endif()
  message(STATUS "Found PROGRESS: ${PROGRESS_LIBRARIES} ${BML_LIBRARIES}")
endif()  

option(DBCSR_OPT "Whether to use DBCSR" OFF)

option(GPUOPT "Whether to use the GPU" OFF)
if(GPUOPT)
  find_package(CUDA REQUIRED)
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_CURRENT_SOURCE_DIR}/../src/Doxyfile.in Doxyfile)
  add_custom_target(doc
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_BINARY_DIR}/Doxyfile)
endif()

set(MYMODULES  myprecision.f90 constants_mod.f90 setuparray.f90 coulombarray.f90 
diagarray.f90 fermicommon.f90 mdarray.f90 neblistarray.f90 
nonoarray.f90 ppotarray.f90 purearray.f90 relaxcommon.f90 restartarray.f90 
sparsearray.f90 spinarray.f90 univarray.f90 virialarray.f90 xboarray.f90 
kspacearray.f90 timer_mod.f90 subgraph.f90 matrixio.f90 sparsemath.f90 
sparsesp2.f90 subgraphsp2.f90 homolumo.f90 openfiles_mod.f90 kernelparser_mod.f90 
latteparser_latte_mod.f90 dmarray.f90 mixer_mod.f90 sp2progress.f90 genXprogress.f90  
bodirectprogress.f90 constraints_mod.f90 latte_lib.f90 wignerarray.f90)

set(DBCSRMODS dbcsr_var_mod.f90 init_dbcsr.f90 shutdown_dbcsr.f90)

set(MYSUBS allocatediag.f90 allocatenono.f90 allocatepure.f90 bodirect.f90 
deallocatediag.f90 deallocatepure.f90 deorthomyrho.f90 deallocatenono.f90 
allocatedm.f90 deallocatedm.f90 initiatedm.f90 adddftbu.f90 xbodm.f90 
hubbardforce.f90 hubfcontr_sp.f90 Invert.f90 
diagmyh.f90 diagmyhprogress.f90 entropy.f90 fermiexpans.f90 genX.f90  latte.F90 orthomyH.f90 
 orthomyrho.f90  pulay.f90  solvematcg.f90 solvematlapack.f90 
sp2T.f90 pulay_spprogress.f90 
sp2fermi.f90 sp2fermi_init.f90 sp2pure.f90 sp2pure_sparse.f90 kdiagmyh.f90 
 kbodirect.f90 sp2gap.f90 sp2gap_setup.f90 readrestart.f90 initrng.f90 
parawrite.f90 parafileopen.f90 wrtcfgs.f90 wrtrestart.f90 pbc.f90 
addqdep.f90 allocatecoulomb.f90 allocatenebarrays.f90  allocatexbo.f90  
atomcharge.f90 avepress.f90 avetemp.f90 bldnewH.f90 bldspinH.f90 
conjgradient.f90 coulomb_ewald.f90 coulomb_oldskool.f90 coulomb_rspace.f90 
coultailcoef.f90 deallocateall.f90 deallocatecoulomb.f90 
deallocatenebarrays.f90  deallocatexbo.f90 
fermiallocate.f90 fermideallocate.f90 
   gershgorin.f90 getbndfil.f90 
getcoule.f90 getdensity.f90 getdeltaq.f90 getdeltaspin.f90 gethdim.f90 
getke.f90 getmaxf.f90 getmdf.f90 getpressure.f90 getrho.f90 getspinE.f90  
 ifrestart.f90 initcoulomb.f90 initialv.f90 initshockcomp.f90 msrelax.f90 
neblists.f90 noelec.f90 norms.f90 nptrescale.f90  nvtrescale.f90 
pairpot.f90 panic.f90 plot_ppot.f90 plot_univ.f90 printsparse.f90 
 propchempot_xbo.f90 propspins_xbo.f90 qconsistency.f90 readcontrols.f90 
readcr.f90 readmdcontroller.f90 readppot.f90  readtb.f90 resetprodhd.f90
rhozero.f90 shockcomp.f90  spinrhodirect.f90 spinrhodirectprogress.f90 stdescent.f90 
summary.f90 tbmd.f90 toteng.f90 univscaling.f90 univscaling_function.f90 univtailcoef.f90  
vdwtailcoef.f90 velverlet.f90  xbo.f90 am.f90 bm.f90  
 wignerd.f90 dfdx.f90 dfda.f90 dfdb.f90 dfdr.f90 dfda_prec.f90 dfdb_prec.f90 dfdr_prec.f90 
 slmmp.f90 tlmmp.f90 dslmmpda.f90 dslmmpdb.f90 dtlmmpda.f90 dtlmmpdb.f90 dwignerddb.f90  
 gradH.f90 kbldnewh.f90  kgetrho.f90 kgradH.f90  
  dunivscaling.f90 kgetdos.f90 
qneutral.f90 shiftH.f90 getmatindlist.f90 nvtlangevin.f90 nvtandersen.f90 
nvtNH.f90 sp2pure_sparse_parallel.f90 sp2pure_sparse_parallel_simple.f90 
sp2pure_subgraph_parallel.f90 gaussrn.f90 nnz.f90 allocatesubgraph.f90 
deallocatesubgraph.f90 dosfit.f90 mofit.f90 getrespf.f90 gendiag.f90 
mofit_plato.f90 pairpot_noneb.f90  
readppottab.f90 pairpottab.f90 hugrescale.f90 avesforhug.f90 getdipole.f90 
fittingoutput.f90 tabtest.f90 kgenX.f90 korthomyH.f90 kdeorthomyrho.f90 
 getforce.f90 masses2symbols.f90 setuptbmd.f90 
wrtrestartlib.f90 deorthomyrhoprogress.f90 orthomyHprogress.f90 readrestartlib.f90 
readppotspline.f90 pairpotspline.f90 errors.f90 initcoulombklist.f90 
  dosfittab.f90 
readppotD.f90 pairpot_plusD.f90  tbfit.f90 fire.f90 angfactor.f90 
build_integral_map.f90 tb_forces.f90 dfdx_HS.f90 genorbitallist.f90 
ktb_forces.f90 genHonsite.f90 gencutofflist.f90 getcorrs.f90 
Canon_DM_PRT.f90 can_resp.f90 getdeltaq_resp.f90  bldnewHS_sp.f90) 

foreach(_m MYMODULES MYSUBS DBCSR)
  set(${_m}_SOURCES)
  foreach(_file ${${_m}})
    list(APPEND ${_m}_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/../src/${_file})
  endforeach()
endforeach()

if(NOT DBCSR_OPT)
  set(DBCSR_SOURCES)
endif()

add_library(latte ${MYMODULES_SOURCES} ${MYSUBS_SOURCES}
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/latte_lib.f90
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/latte_c_bind.f90
)
target_compile_options(latte PRIVATE -cpp)
install(TARGETS latte EXPORT LATTE_Targets LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR} ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(EXPORT LATTE_Targets FILE LATTE_Targets.cmake NAMESPACE LATTE:: DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTE)
include(CMakePackageConfigHelpers)
configure_file(LATTEConfig.cmakein ${CMAKE_CURRENT_BINARY_DIR}/LATTEConfig.cmake @ONLY)
write_basic_package_version_file("LATTEConfigVersion.cmake" VERSION ${PROJECT_VERSION} COMPATIBILITY ExactVersion)
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/LATTEConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/LATTEConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/LATTE)
add_library(LATTE::latte ALIAS latte)
target_link_libraries(latte PRIVATE
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
  target_link_libraries(latte PUBLIC OpenMP::OpenMP_Fortran)
endif()
target_include_directories(latte PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${PROGRESS_INCLUDE_DIRS} ${BML_INCLUDE_DIRS})

#only "ON" and "OFF" support
foreach(_def GPUOPT DBCSR_OPT DO_MPI PROGRESS)
  if(${${_def}})
    set(${_def} ON)
  else()
    set(${_def} OFF)
  endif()
endforeach()

target_compile_definitions(latte PRIVATE XSYEV ${PRECISION}PREC GPU${GPUOPT} DBCSR_${DBCSR_OPT} MPI_${DO_MPI} PROGRESS${PROGRESS})

if(DO_MPI)
  set(LATTE LATTE_MPI_${PRECISION})
else()
  set(LATTE LATTE_${PRECISION})
endif()
add_executable(${LATTE} ${CMAKE_CURRENT_SOURCE_DIR}/../src/latte.F90)
target_link_libraries(${LATTE} latte ${PROGRESS_LIBRARIES})
target_include_directories(${LATTE} PRIVATE ${MPI_Fortran_INCLUDE_PATH} ${PROGRESS_INCLUDE_DIRS} ${BML_INCLUDE_DIRS})
target_compile_definitions(${LATTE} PRIVATE XSYEV ${PRECISION}PREC GPU${GPUOPT} DBCSR_${DBCSR_OPT} MPI_${DO_MPI} PROGRESS${PROGRESS})
install(TARGETS ${LATTE} DESTINATION ${CMAKE_INSTALL_BINDIR})
