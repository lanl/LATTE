#
# LATTE Makefile
#

include makefile.CHOICES

PROGRAMS = LATTE_SINGLE LATTE_DOUBLE LATTEGPU_SINGLE LATTEGPU_DOUBLE LATTE_DBCSR_DOUBLE LATTE_DBCSR_SINGLE

all : 
ifeq ($(GPUOPT),ON)
	(cd MATRIX; make; cd ..)
endif
	(rm liblatte.a; cd src; make; cd ..)

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
