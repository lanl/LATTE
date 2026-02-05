source vars
MY_PATH=`pwd`
echo $MY_PATH
cd $MY_PATH/lammps/lib/latte
rm includelink; ln -s $MY_PATH/LATTE/src includelink
rm liblink; ln -s $MY_PATH/LATTE liblink
rm filelink.o; ln -s $MY_PATH/LATTE/src/latte_c_bind.o filelink.o
cd $MY_PATH
MAGMA_PATH=/projects/shared/spack/opt/spack/linux-ubuntu22.04-zen4/gcc-11.4.0/magma-2.9.0-brkm5qdkhezbomtt36p54a7iyagrbcso
METIS_PATH=/projects/shared/spack/opt/spack/linux-ubuntu22.04-zen4/gcc-11.4.0/metis-5.1.0-pdxbuygwakvh2qaxiskhf3qkyyk5cmiv
#Construct Makefile for lammps

echo "progress_PATH = ${MY_PATH}/qmd-progress" > Makefile.lammps
echo "bml_PATH = ${MY_PATH}/bml" >> Makefile.lammps
echo "latte_SYSLIB += -L${CUDA_PATH}/lib64/ -lcublas -lcudart -lcusolver" >> Makefile.lammps
echo "latte_SYSLIB += -I${CUDA_PATH}/include" >> Makefile.lammps
echo "latte_SYSLIB += ../../lib/latte/filelink.o -llatte" >> Makefile.lammps
echo "latte_SYSLIB += -L${MAGMA_PATH}/lib -lmagma -lm -lgfortran" >> Makefile.lammps 
echo "latte_SYSLIB += -L${METIS_PATH}/lib -lmetis" >> Makefile.lammps 
echo "latte_SYSINC  +=  -I${MY_PATH}/bml/install/include -I${MY_PATH}/qmd-progress/install/include" >> Makefile.lammps 
echo "latte_SYSINC  += -Wl,--copy-dt-needed-entries" >> Makefile.lammps 
echo "latte_SYSINC += -I${METIS_PATH}/include " >> Makefile.lammps 
echo "latte_SYSLIB  += -L${MY_PATH}/bml/install/lib -lprogress -L${MY_PATH}/qmd-progress/install/lib -lbml_fortran -lbml" >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp " >> Makefile.lammps
#echo "latte_SYSLIB += -L${MKLROOT}/lib -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64" >> Makefile.lammps
#echo "latte_SYSLIB += -fopenmp -lpthread -liomp5" >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp -lpthread " >> Makefile.lammps
echo "latte_SYSLIB += -Wl,--copy-dt-needed-entries " >> Makefile.lammps
echo "latte_SYSLIB += -llapack -lblas " >> Makefile.lammps

cp Makefile.lammps ./lammps/lib/latte/
cd ./lammps/src/ ; make yes-molecule; make yes-latte; make -j serial
cd $MY_PATH
