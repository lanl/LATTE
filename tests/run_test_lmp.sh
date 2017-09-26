#!/bin/bash

# Script to test the LAMMPS-LATTE interface program.

set -e                                          # This will exit the script if there is any error
MY_PATH=`pwd`                                   # Capturing the local path of the folder where we are running.

RUN=$HOME"/lammps/src/lmp_serial"  #EXAALT version of Lammps program

for name in 0scf 2scf fullscf ; do

  INLATTEFILE="latte."$name".in"
  INLAMMPSFILE="in."$name
  REF="energy."$name".out"
  COORDS="data."$name".lmp"

  cp  ./tests/tests_lmp/$INLATTEFILE latte.in
  cp  ./tests/tests_lmp/$INLAMMPSFILE .
  cp  ./tests/tests_lmp/$REF .
  cp  ./tests/tests_lmp/$COORDS .

  echo -e "\nTesting for "$name" \n"

  time $RUN < $INLAMMPSFILE >  out
  grep Energy out | sed -e s/"Total Energy ="/""/g >  energy.out
  echo ""

  grep Energy out | sed -e s/"PAR"/$STRR/g  >  input_tmp.in
  python ./tests/test-energy.py --reference $REF --current energy.out --reltol 0.000001

done

# Tests for geometry optimizations 

for name in opt ; do

  INLATTEFILE="latte."$name".in"
  INLAMMPSFILE="in."$name
  REF="energy."$name".out"
  COORDS="data."$name".lmp"

  cp  ./tests/tests_lmp/$INLATTEFILE latte.in
  cp  ./tests/tests_lmp/$INLAMMPSFILE .
  cp  ./tests/tests_lmp/$REF .
  cp  ./tests/tests_lmp/$COORDS .

  echo -e "\nTesting for "$name" \n"

  time $RUN < $INLAMMPSFILE >  out
  grep -A 16 TotEng out | sed -e s/"TotEng"/"0.0"/g >  energy.out
  echo ""
  
  python ./tests/test-energy.py --reference $REF --current energy.out --reltol 0.0000001

done

echo -e "\nEnd of run and test"
