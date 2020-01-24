#!/bin/bash

# Script to test the LAMMPS-LATTE interface program.

cp latte.in latte.tmp

set -e                                          # This will exit the script if there is any error
MY_PATH=`pwd`                                   # Capturing the local path of the folder where we are running.

RUN=$HOME"/software/lammps/src/lmp_serial"  #EXAALT version of Lammps program

for name in 0scf 2scf fullscf 0scf.wrtrestart 0scf.rdrestart ; do

  INLATTEFILE="latte."$name".in"
  INLAMMPSFILE="in."$name
  REF="energy."$name".out"
  COORDS="data."$name".lmp"

  if [ "$name" == "0scf.rdrestart" ] ; then  
	cp restart.latte.10.dat restart.latte.dat	
  fi

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

for name in opt opt.boxrel ; do

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

cp latte.tmp latte.in
rm *.in *.out out *.lmp in.* log.* restart.* latte.tmp

echo -e "\nEnd of run and test"
