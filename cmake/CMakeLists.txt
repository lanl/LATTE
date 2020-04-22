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

set(MYMODULES myprecision.F90 constants_mod.F90 setuparray.F90 coulombarray.F90
diagarray.F90 fermicommon.F90 mdarray.F90 neblistarray.F90
nonoarray.F90 ppotarray.F90 purearray.F90 relaxcommon.F90 restartarray.F90
sparsearray.F90 spinarray.F90 univarray.F90 virialarray.F90 xboarray.F90
kspacearray.F90 timer_mod.F90 subgraph.F90 matrixio.F90 sparsemath.F90
sparsesp2.F90 subgraphsp2.F90 homolumo.F90 openfiles_mod.F90 kernelparser_mod.F90
latteparser_latte_mod.F90 mixer_mod.F90 sp2progress.F90 genXprogress.F90
bodirectprogress.F90 constraints_mod.F90 latte_lib.F90)

set(DBCSRMODS dbcsr_var_mod.F90 init_dbcsr.F90 shutdown_dbcsr.F90)

set(MYSUBS allocatediag.F90 allocatenono.F90 allocatepure.F90 bodirect.F90
deallocatediag.F90 deallocatepure.F90 deorthomyrho.F90 deallocatenono.F90
diagmyh.F90 entropy.F90 fermiexpans.F90 genX.F90  latte.F90 orthomyH.F90
 orthomyrho.F90  pulay.F90 pulay_sp.F90 solvematcg.F90 solvematlapack.F90
sp2T.F90 pulay_spprogress.F90
sp2fermi.F90 sp2fermi_init.F90 sp2pure.F90 sp2pure_sparse.F90 kdiagmyh.F90
 kbodirect.F90 sp2gap.F90 sp2gap_setup.F90 readrestart.F90 initrng.F90
parawrite.F90 parafileopen.F90 wrtcfgs.F90 wrtrestart.F90 pbc.F90
addqdep.F90 allocatecoulomb.F90 allocatenebarrays.F90  allocatexbo.F90 
atomcharge.F90 avepress.F90 avetemp.F90 bldnewH.F90 bldspinH.F90
conjgradient.F90 coulomb_ewald.F90 coulomb_oldskool.F90 coulomb_rspace.F90
coultailcoef.F90 deallocateall.F90 deallocatecoulomb.F90
deallocatenebarrays.F90  deallocatexbo.F90 fcoulnono.F90 fcoulnono_sp.F90
fermiallocate.F90 fermideallocate.F90
fspinnono.F90 fspinnono_sp.F90  gershgorin.F90 getbndfil.F90
getcoule.F90 getdensity.F90 getdeltaq.F90 getdeltaspin.F90 gethdim.F90
getke.F90 getmaxf.F90 getmdf.F90 getpressure.F90 getrho.F90 getspinE.F90 
 ifrestart.F90 initcoulomb.F90 initialv.F90 initshockcomp.F90 msrelax.F90
neblists.F90 noelec.F90 norms.F90 nptrescale.F90  nvtrescale.F90
pairpot.F90 panic.F90 plot_ppot.F90 plot_univ.F90 printsparse.F90
 propchempot_xbo.F90 propspins_xbo.F90 qconsistency.F90 readcontrols.F90
readcr.F90 readmdcontroller.F90 readppot.F90  readtb.F90 resetprodhd.F90
rhozero.F90 shockcomp.F90  spinrhodirect.F90 stdescent.F90
summary.F90 tbmd.F90 toteng.F90 univscaling_function.F90 univtailcoef.F90 
vdwtailcoef.F90 velverlet.F90  xbo.F90 am.F90 bm.F90 factorial.F90
 slmmp.F90 tlmmp.F90 wignerd.F90 dfda.F90 dfdb.F90 dfdr.F90
dunivscaling_function.F90 dslmmpda.F90 dslmmpdb.F90 dtlmmpda.F90
dtlmmpdb.F90 dwignerddb.F90 gradH.F90 kbldnewh.F90  kgetrho.F90 kgradH.F90 
bldnewHS_sp.F90 gradH_sp.F90 univscaling.F90 dunivscaling.F90 kgetdos.F90
qneutral.F90 shiftH.F90 getmatindlist.F90 nvtlangevin.F90 nvtandersen.F90
nvtNH.F90 sp2pure_sparse_parallel.F90 sp2pure_sparse_parallel_simple.F90
sp2pure_subgraph_parallel.F90 gaussrn.F90 nnz.F90 allocatesubgraph.F90
deallocatesubgraph.F90 dosfit.F90 mofit.F90 getrespf.F90 gendiag.F90
mofit_plato.F90 ppfit.F90 pairpot_noneb.F90 allfit.F90
readppottab.F90 pairpottab.F90 hugrescale.F90 avesforhug.F90 getdipole.F90
fittingoutput.F90 tabtest.F90 kgenX.F90 korthomyH.F90 kdeorthomyrho.F90
kfcoulnono.F90 kpulay.F90 getforce.F90 masses2symbols.F90 setuptbmd.F90
wrtrestartlib.F90 deorthomyrhoprogress.F90 orthomyHprogress.F90 readrestartlib.F90
readppotspline.F90 pairpotspline.F90 errors.F90 initcoulombklist.F90
flcnnono.F90 flcnnono_sp.F90 kflcnnono.F90 dosfittab.F90 
readppotD.F90 pairpot_plusD.F90)

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
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/latte_lib.F90
  ${CMAKE_CURRENT_SOURCE_DIR}/../src/latte_c_bind.F90
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
