#
# LATTE Makefile
#

include makefile.CHOICES

PROGRAMS = LATTE_SINGLE LATTE_DOUBLE LATTEGPU_SINGLE LATTEGPU_DOUBLE LATTE_DBCSR_DOUBLE LATTE_DBCSR_SINGLE

MY_PATH=$(shell pwd)

all : 
ifeq ($(GPUOPT),ON)
	(cd MATRIX; make; cd ..)
endif
	(rm liblatte.a; cd src; make; cd ..)

lammps : 
	(rm liblatte.a; cd src; make; cd ..)
	(cd $(HOME)/lammps/src; touch fix_latte.cpp; make serial; cd -)
	
src : 
	(rm liblatte.a; cd src; make; cd ..)

docs :
	(cd ./src; doxygen Doxyfile.in)
	(cd ./doc/latex/ ; make ; cd ../../)
	(cp ./doc/latex/refman.pdf ./Manual/)

test : 
	(./tests/run_test.sh)

test_lmp : 
	(./tests/run_test_lmp.sh)

cov : 
	(./tests/run_cov)

matrix : 
ifeq ($(GPUOPT),ON)
	(cd MATRIX; make; cd ..)
endif

clean : 
	(cd MATRIX; make clean; cd ..)
	(cd src; make clean; cd ..)
	rm *.a

veryclean : 
	(cd MATRIX; make clean; cd ..)
	(cd src; make clean; cd ..)
	rm $(PROGRAMS) *.a
