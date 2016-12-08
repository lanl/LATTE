include ../makefile.CHOICES

DEBUG=
PROF=
CUDA_OPT_LEVEL=3
CC_OPT_LEVEL=0
INCLUDES = -I. -IKernels

# Set number size based on precision
ifeq ($(PRECISION),SINGLE)
  REALSIZE = 4
else
  REALSIZE = 8
endif

CUDA_CC= nvcc -arch $(GPU_ARCH) -O${CUDA_OPT_LEVEL} ${PROF} ${DEBUG} ${INCLUDES} -DREALSIZE=$(REALSIZE)

KERNEL_SOURCES=$(wildcard Kernels/*.cu)
KERNEL_OBJECTS=$(KERNEL_SOURCES:.cu=.$(REALSIZE).o)

#MATRIX_SOURCES=$(wildcard M_*.cu) $(wildcard F_*.cu) SP2Pure.cu SP2Fermi.cu SolveMatrixCG.cu Allocate.cu TestMultiply.cu
MATRIX_SOURCES=$(wildcard M_*.cu) $(wildcard F_*.cu) SP2Pure.cu SP2Fermi.cu SolveMatrixCG.cu Allocate.cu runmatmult.cu
MATRIX_CUDA_OBJECTS=$(MATRIX_SOURCES:.cu=.cuda.$(REALSIZE).o)

all: libmatrix_cuda.$(PRECISION).a

libmatrix_cuda.$(PRECISION).a: $(MATRIX_CUDA_OBJECTS) $(KERNEL_OBJECTS)
	ar r libmatrix_cuda.$(PRECISION).a $(MATRIX_CUDA_OBJECTS) $(KERNEL_OBJECTS)

%.cuda.$(REALSIZE).o: %.cu
	${CUDA_CC} -c $< -o $@ 

%.$(REALSIZE).o: %.cu
	${CUDA_CC} -c $< -o $@

clean:
	- rm *~ *.o *.a Kernels/*~ Kernels/*.o
