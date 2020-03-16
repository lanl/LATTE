#!/bin/bash

#SBATCH -n36 -N1 #Number of cores and nodes
##SBATCH --ntasks-per-node=36
#SBATCH --time=10:00:00
#SBATCH --partition=scaling   
## OMP_NUM_THREADS controls the number of threads your application use
## This variable cat be set by the following command :


## Run command
#./run.sh | tee out 

export OMPI_MCA_opal_paffinity_alone=0
#export MKL_DISABLE_FAST_MM=1
export OMP_NUM_THREADS=36
export NSLOTS=200

module purge
module load gcc/6.4.0
#module load openmpi
module load mkl
module load cmake

#$HOME/ecp/lammps/src/lmp_serial < run.lammps >  out 

cp latte.in-kernel latte.in
sed -i "s/NORECS= 0/NORECS= 2/g" latte.in

../../LATTE-ecp/LATTE_DOUBLE < latte.in  > out

#for i in 1 2 4 8 12
##for i in 4
##do
##  echo $i
##  cp latte.in-kernel latte.in
##  sed -i "s/NORECS= 0/NORECS= $i/g" latte.in
##
##  ../../LATTE-ecp/LATTE_DOUBLE < latte.in  > out
##  grep 'Data ' out > energy.kernel$i
##  grep 'MDITER ' out > res.kernel$i
##done


##
###echo "done"
##
##mkdir dt1.0
##mv energy.* dt1.0
##mv res.* dt1.0
##
##
##for i in 1 2 4 8 12
##do
##  echo $i
##  cp latte.in-kernel-dt latte.in
##  sed -i "s/NORECS= 0/NORECS= $i/g" latte.in
##
##  ../../LATTE-ecp/LATTE_DOUBLE < latte.in  > out
##  grep 'Data ' out > energy.kernel$i
##  grep 'MDITER ' out > res.kernel$i
##done
