#
# LATTE Makefile
#

include makefile.CHOICES

PROGRAMS = LATTE_SINGLE LATTE_DOUBLE LATTEGPU_SINGLE LATTEGPU_DOUBLE LATTE_DBCSR_DOUBLE LATTE_DBCSR_SINGLE

all : 
ifeq ($(GPUOPT),ON)
	(cd MATRIX; make; cd ..)
endif
	(cd src; make; cd ..)

src : 
	(cd src; make; cd ..)

docs :
	(doxygen ./src/Doxyfile.in)

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
