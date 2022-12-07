#!/bin/sh

topdir=$PWD
echo $topdir
if [ !-d "build" ]; then
  mkdir build
fi

cd build

# debug
#cmake ../cmake -DCMAKE_INSTALL_PREFIX=$topdir/install -DCMAKE_Fortran_FLAGS='-O0 -g -traceback -check all -debug all -check all -fpe0' #-DOPENMP=${CMAKE_WITH_PROGRESS} 

# release
cmake ../cmake -DCMAKE_INSTALL_PREFIX=$topdir/install -DCMAKE_Fortran_FLAGS='-O3 -fpp' #-DOPENMP=${CMAKE_WITH_PROGRESS} 

VERBOSE=2 make
make install



#cd ~/path/to/lammps_compute_pace/build_latte/
#rm CMakeFiles/lammps.dir/vast/home/zhy/ecp/lammps/src/LATTEQEQ/fix_latteqeq.cpp.o
#~/path/to/lammps_compute_pace//build_pace_latteqeq.sh 
#make -j

cd ~/ecp/mylammps/build_latte
rm CMakeFiles/lammps.dir/vast/home/zhy/ecp/lammps/src/LATTEQEQ/fix_latteqeq.cpp.o
~/ecp/mylammps/build_pace_latteqeq.sh 


