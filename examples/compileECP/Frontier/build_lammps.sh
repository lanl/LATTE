MY_PATH=`pwd`
echo $MY_PATH
cd $MY_PATH/lammps/lib/latte
rm includelink; ln -s $MY_PATH/LATTE/src includelink
rm liblink; ln -s $MY_PATH/LATTE liblink
rm filelink.o; ln -s $MY_PATH/LATTE/src/latte_c_bind.o filelink.o
cd $MY_PATH

#Construct Makefile for lammps

echo "progress_PATH = ${MY_PATH}/qmd-progress" > Makefile.lammps
echo "bml_PATH = ${MY_PATH}/bml" >> Makefile.lammps
echo "latte_SYSLIB += ../../lib/latte/filelink.o -llatte" >> Makefile.lammps
echo "latte_SYSINC  +=  -I${MY_PATH}/bml/install/include -I${MY_PATH}/qmd-progress/install/include" >> Makefile.lammps 
echo "latte_SYSLIB  += -L${MY_PATH}/bml/install/lib64 -lprogress -L${MY_PATH}/qmd-progress/install/lib64 -lbml_fortran -lbml" >> Makefile.lammps
echo "latte_SYSINC  +=  -I${OLCF_MAGMA_ROOT}/include -I${ROCM_PATH}/include" >> Makefile.lammps
echo "latte_SYSLIB  += -L${BML_PATH}/install/lib64 -lprogress -L${PROGRESS_PATH}/install/lib64 -lbml_fortran -lbml" >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp -L${OLCF_MAGMA_ROOT}/lib  -lmagma -L${ROCM_PATH}/lib -lrocblas -lamdhip64 -lrocsolver" >> Makefile.lammps
echo "latte_SYSLIB += -lm -lgfortran" >> Makefile.lammps
echo "latte_SYSLIB += -fopenmp -lpthread -L${OLCF_OPENBLAS_ROOT}/lib -L${OLCF_OPENBLAS_ROOT}/lib  -lopenblas" >> Makefile.lammps
echo "latte_SYSINC += -g -I/opt/rocm-4.5.2/include -O2 -DNDEBUG -fopenmp -fopenmp" >> Makefile.lammps
echo "latte_SYSINC += -I/opt/rocm-5.2.0/include -I${PROGRESS_PATH}/build/src -I${BML_PATH}/install/include" >> Makefile.lammps
echo "latte_SYSLIB += -g -ffree-form -ffree-line-length-none -ffree-line-length-none -ffree-line-length-none -O2 -g -ffree-line-length-none -fopenmp -fopenmp" >> Makefile.lammps

cp Makefile.lammps ./lammps/lib/latte/
cd ./lammps/src/ ; make yes-molecule; make yes-latte; make -j serial
cd $MY_PATH
