source vars
MY_PATH=`pwd`
echo $MY_PATH
cd $MY_PATH/lammps/lib/latte
rm includelink; ln -s $MY_PATH/LATTE/src includelink
rm liblink; ln -s $MY_PATH/LATTE liblink
rm filelink.o; ln -s $MY_PATH/LATTE/src/latte_c_bind.o filelink.o
cd $MY_PATH

LBML_PATH=${BML_PATH}
LMAGMA_PATH=${MAGMA_PATH}
LMETIS_PATH=${METIS_PATH}
LOPENBLAS_PATH=${OPENBLAS_PATH}
LCUDA_PATH=${CUDA_PATH}
LPROGRESS_PATH=${PROGRESS_PATH}

#Construct Makefile for lammps

echo "progress_PATH = ${MY_PATH}/qmd-progress" > Makefile.lammps
echo "bml_PATH = ${LBML_PATH}/" >> Makefile.lammps
echo "latte_SYSLIB += -L${LCUDA_PATH}/lib64/ -lcublas -lcudart -lcusolver" >> Makefile.lammps
echo "latte_SYSLIB += -I${LCUDA_PATH}/include" >> Makefile.lammps
echo "latte_SYSLIB += ../../lib/latte/filelink.o -llatte" >> Makefile.lammps
echo "latte_SYSLIB += -L${LMAGMA_PATH}/lib -lmagma -lm -lgfortran" >> Makefile.lammps 
#echo "latte_SYSLIB += -L${LMETIS_PATH}/lib -lmetis" >> Makefile.lammps 
echo "latte_SYSINC  +=  -I${LBML_PATH}/include -I${LPROGRESS_PATH}/include" >> Makefile.lammps 
echo "latte_SYSINC  += -Wl,--copy-dt-needed-entries" >> Makefile.lammps 
echo "latte_SYSINC += -I${LMETIS_PATH}/include " >> Makefile.lammps 
echo "latte_SYSINC += -I${LOPENBLAS_PATH}/include " >> Makefile.lammps 
echo "latte_SYSLIB += -L${LOPENBLAS_PATH}/lib -lopenblas " >> Makefile.lammps 
echo "latte_SYSLIB  += -L${LBML_PATH}/lib -lprogress -L${LPROGRESS_PATH}/lib -lbml_fortran -lbml" >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp " >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp -lpthread " >> Makefile.lammps
echo "latte_SYSLIB += -Wl,--copy-dt-needed-entries " >> Makefile.lammps
echo "latte_SYSLIB += -llapack -lblas " >> Makefile.lammps

cp Makefile.lammps ./lammps/lib/latte/
cd ./lammps/src/ ; make yes-molecule; make yes-latte; make -j serial
cd $MY_PATH
