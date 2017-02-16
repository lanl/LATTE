#!/bin/bash

# Script to test the LAMMPS-LATTE interface program.

set -e                                          # This will exit the script if there is any error
MY_PATH=`pwd`                                   # Capturing the local path of the folder where we are running.

RUN=$HOME"/lammps/src/lmp_serial"  #EXAALT version of Lammps program  

for name in 0scf fullscf ; do

  INLATTEFILE="latte."$name".in"
  INLAMMPSFILE="in."$name
  REF="energy."$name".out"
  COORDS="data."$name".lmp"

  cp  ./example_lmp/tests/$INLATTEFILE latte.in
  cp  ./example_lmp/tests/$INLAMMPSFILE .
  cp  ./example_lmp/tests/$REF .
  cp  ./example_lmp/tests/$COORDS .

  echo -e "\nTesting for "$name" \n"

  time $RUN < $INLAMMPSFILE >  out
  grep Energy out | sed -e s/"Total Energy ="/""/g >  energy.out
  echo ""

  grep Energy out | sed -e s/"PAR"/$STRR/g  >  input_tmp.in
  python ./example_lmp/test-energy.py --reference $REF --current energy.out --reltol 0.0000001

done

echo -e "\nEnd of run and test"
