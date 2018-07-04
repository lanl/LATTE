include ../makefile.CHOICES

PPFLAGS = -DXSYEVD -D$(PRECISION)PREC -DGPU$(GPUOPT) -DDBCSR_$(DBCSR_OPT) -DMPI_$(MPIOPT)

all:
# 	$(FC) $(FFLAGS) $(PPFLAGS) timer_cpu_time.f90 -o timer_cpu_time
	$(FC) timer_cpu_time.f90 -o timer_cpu_time
