#!/bin/bash

# Script to test the LATTE program.

MY_PATH=`pwd`                                   # Capturing the local path of the folder where we are running.

RUN="./LATTE_DOUBLE"  # LATTE program executable

echo -e "\nTesting LATTE with new (latte.in) input files \n"

mv latte.in latte.in.tmp

set -e                                          # This will exit the script if there is any error


# Testing for single point calculations:

for name in single.point single.point.noelec single.point.rspace ; do  

  INLATTEFILE="latte."$name".in"
  REF="energy."$name".out"
  COORDS=$name".dat"

  cp  ./tests/$INLATTEFILE latte.in
  cp  ./tests/$REF .
  cp  ./tests/$COORDS ./bl/inputblock.dat

  echo -e "\nTesting for "$name" \n"

  time $RUN > out
  ENERG=`grep "FREE" out | awk 'NF>1{print $5}'`
  echo $ENERG > energy.out

  python ./tests/test-energy.py --reference $REF --current energy.out --reltol 0.00001

done 

# Testing geometry optimizations:

for name in opt opt.cg ; do

  INLATTEFILE="latte."$name".in"
  REF="monitorrelax."$name".xyz"
  COORDS=$name".dat"

  cp  ./tests/$INLATTEFILE latte.in
  cp  ./tests/$REF .
  cp  ./tests/$COORDS ./bl/inputblock.dat

  echo -e "\nTesting for "$name" \n"

  time $RUN > out

  python ./tests/test-optim.py --reference $REF --current monitorrelax.xyz --reltol 0.00001

done

# Testing for MD simulations:

for name in 0scf 2scf fullscf fullscf.etemp sp2 sp2.sparse fullscf.nvt \
       	fullscf.npt fullscf.vdw fullscf.spin fullscf.kon ; do

  INLATTEFILE="latte."$name".in"
  REF="energy."$name".out"
  COORDS=$name".dat"

  cp  ./tests/$INLATTEFILE latte.in
  cp  ./tests/$REF .
  cp  ./tests/$COORDS ./bl/inputblock.dat

  echo -e "\nTesting for "$name" \n"

  time $RUN > out
  grep "Data" out | sed 's/Data/ /g' | awk 'NF>1{print $2}' > energy.out
  echo ""

  python ./tests/test-energy.py --reference $REF --current energy.out --reltol 0.00001

done


# Testing with the usual latte input method:

rm latte.in

echo -e "\nTesting LATTE with original input files \n"

 for name in 0scf fullscf sp2 ; do

   CONTROL="control."$name".in"
   MDCONTROLLER="MDcontroller."$name
   REF="energy."$name".out"
   COORDS=$name".dat"

   cp  ./tests/$CONTROL ./TBparam/control.in
   cp  ./tests/$MDCONTROLLER MDcontroller
   cp  ./tests/$REF .
   cp  ./tests/$COORDS ./bl/inputblock.dat

   echo -e "\nTesting for "$name" \n"

   time $RUN > out
   grep "Data" out | sed 's/Data/ /g' | awk 'NF>1{print $2}' > energy.out
   echo ""

   grep Energy out | sed -e s/"PAR"/$STRR/g  >  input_tmp.in
   python ./tests/test-energy.py --reference $REF --current energy.out --reltol 0.00001
 
 done

mv latte.in.tmp latte.in
rm out 



echo -e "\nEnd of run and test"
