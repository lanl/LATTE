#
# This is the main makefile for the serial version of LATTE
#
.SUFFIXES: .fortran90 .f90

include ../makefile.CHOICES

PPFLAGS = -DXSYEV -D$(PRECISION)PREC -DGPU$(GPUOPT) -DDBCSR_$(DBCSR_OPT) -DMPI_$(MPIOPT)
PPFLAGS += -DPROGRESS$(PROGRESS)

# Flags for linking GPU libraries
GPU_LINKFLAG = -L../MATRIX -lgfortran -lstdc++
GPU_LIB = -L../MATRIX -lmatrix_cuda.$(PRECISION) $(GPU_CUDA_LIB)


MYMODULES = myprecision.f90 constants_mod.f90 setuparray.f90 coulombarray.f90 \
diagarray.f90 fermicommon.f90 mdarray.f90 neblistarray.f90 \
nonoarray.f90 ppotarray.f90 purearray.f90 relaxcommon.f90 restartarray.f90 \
sparsearray.f90 spinarray.f90 univarray.f90 virialarray.f90 xboarray.f90 \
kspacearray.f90 timer_mod.f90 subgraph.f90 matrixio.f90 sparsemath.f90 \
sparsesp2.f90 subgraphsp2.f90 homolumo.f90 openfiles_mod.f90 kernelparser_mod.f90 \
latteparser_latte_mod.f90 mixer_mod.f90 latte_lib.f90

DBCSRMODS = dbcsr_var_mod.f90 init_dbcsr.f90 shutdown_dbcsr.f90

MYSUBS = allocatediag.f90 allocatenono.f90 allocatepure.f90 bodirect.f90 \
deallocatediag.f90 deallocatepure.f90 deorthomyrho.f90 deallocatenono.f90 \
diagmyh.f90 entropy.f90 fermiexpans.f90 genX.f90  latte.f90 orthomyH.f90 \
 orthomyrho.f90  pulay.f90 pulay_sp.f90 solvematcg.f90 solvematlapack.f90 \
sp2T.f90 \
sp2fermi.f90 sp2fermi_init.f90 sp2pure.f90 sp2pure_sparse.f90 kdiagmyh.f90 \
 kbodirect.f90 sp2gap.f90 sp2gap_setup.f90 readrestart.f90 initrng.f90 \
parawrite.f90 parafileopen.f90 wrtcfgs.f90 wrtrestart.f90 pbc.f90 \
addqdep.f90 allocatecoulomb.f90 allocatenebarrays.f90  allocatexbo.f90  \
atomcharge.f90 avepress.f90 avetemp.f90 bldnewH.f90 bldspinH.f90 \
conjgradient.f90 coulomb_ewald.f90 coulomb_oldskool.f90 coulomb_rspace.f90 \
coultailcoef.f90 deallocateall.f90 deallocatecoulomb.f90 \
deallocatenebarrays.f90  deallocatexbo.f90 fcoulnono.f90 fcoulnono_sp.f90 \
fermiallocate.f90 fermideallocate.f90 \
fspinnono.f90 fspinnono_sp.f90  gershgorin.f90 getbndfil.f90 \
getcoule.f90 getdensity.f90 getdeltaq.f90 getdeltaspin.f90 gethdim.f90 \
getke.f90 getmaxf.f90 getmdf.f90 getpressure.f90 getrho.f90 getspinE.f90  \
 ifrestart.f90 initcoulomb.f90 initialv.f90 initshockcomp.f90 msrelax.f90 \
neblists.f90 noelec.f90 norms.f90 nptrescale.f90  nvtrescale.f90 \
pairpot.f90 panic.f90 plot_ppot.f90 plot_univ.f90 printsparse.f90 \
 propchempot_xbo.f90 propspins_xbo.f90 qconsistency.f90 readcontrols.f90 \
readcr.f90 readmdcontroller.f90 readppot.f90  readtb.f90 resetprodhd.f90 \
rhozero.f90 shockcomp.f90  spinrhodirect.f90 stdescent.f90 \
summary.f90 tbmd.f90 toteng.f90 univscaling_function.f90 univtailcoef.f90  \
vdwtailcoef.f90 velverlet.f90  xbo.f90 am.f90 bm.f90 factorial.f90 \
 slmmp.f90 tlmmp.f90 wignerd.f90 dfda.f90 dfdb.f90 dfdr.f90 \
dunivscaling_function.f90 dslmmpda.f90 dslmmpdb.f90 dtlmmpda.f90 \
dtlmmpdb.f90 dwignerddb.f90 gradH.f90 kbldnewh.f90  kgetrho.f90 kgradH.f90  \
bldnewHS_sp.f90 gradH_sp.f90 univscaling.f90 dunivscaling.f90 kgetdos.f90 \
qneutral.f90 shiftH.f90 getmatindlist.f90 nvtlangevin.f90 nvtandersen.f90 \
nvtNH.f90 sp2pure_sparse_parallel.f90 sp2pure_sparse_parallel_simple.f90 \
sp2pure_subgraph_parallel.f90 gaussrn.f90 nnz.f90 allocatesubgraph.f90 \
deallocatesubgraph.f90 dosfit.f90 mofit.f90 getrespf.f90 gendiag.f90 \
mofit_plato.f90 ppfit.f90 pairpot_noneb.f90 allfit.f90 \
readppottab.f90 pairpottab.f90 hugrescale.f90 avesforhug.f90 getdipole.f90 \
fittingoutput.f90 tabtest.f90 kgenX.f90 korthomyH.f90 kdeorthomyrho.f90 \
kfcoulnono.f90 kpulay.f90 getforce.f90 masses2symbols.f90 setuptbmd.f90 \
wrtrestartlib.f90 readrestartlib.f90 

# dependencies
# DEP_FILE = makefile.dep
# include $(DEP_FILE)

XMYMODULES = 		$(MYMODULES:.f90=.o)
XDBCSRMODS = 		$(DBCSRMODS:.f90=.o)
XMYSUBS = 		$(MYSUBS:.f90=.o)


ifeq ($(DBCSR_OPT), ON)

latte:		$(XMYMODULES) $(XDBCSRMODS) $(XMYSUBS)
		$(FCL) -o ../LATTE_DBCSR_$(PRECISION) $(LINKFLAG) $(XMYMODULES) $(XDBCSRMODS) $(XMYSUBS) $(DBCSR_LIB) $(LIB) $(DBCSR_MOD)

else
ifeq ($(GPUOPT), ON)

latte:		$(XMYMODULES) $(XMYSUBS)
		$(FCL) -o ../LATTEGPU_$(PRECISION) $(GPU_LINKFLAG) $(XMYMODULES) $(XMYSUBS) $(LIB) $(GPU_LIB)

else
ifeq ($(MPIOPT), ON)

latte: 		$(XMYMODULES) $(XMYSUBS)
		$(FCL) -o ../LATTE_MPI_$(PRECISION) $(LINKFLAG) $(XMYMODULES)  $(XMYSUBS) $(LIB)

else
latte:		$(XMYMODULES) $(XMYSUBS)
		$(FCL) -o ../LATTE_$(PRECISION) $(LINKFLAG) $(XMYMODULES) $(XMYSUBS) $(LIB)

endif
endif
endif


ifeq ($(MAKELIB), ON)
		$(AR) ../liblatte.a latte_lib.o $(XMYMODULES) $(XMYSUBS)
		$(RANLIB) ../liblatte.a
		$(FC) $(FFLAGS) -c latte_c_bind.f90
endif



ifeq ($(DBCSR_OPT), ON)

$(XMYMODULES): 	        $(MYMODULES)
			$(FC) $(FFLAGS) $(PPFLAGS) -c $(MYMODULES) $(DBCSR_MOD)

$(XDBCSRMODS): 	        $(DBCSRMODS)
			$(FC) $(FFLAGS) $(PPFLAGS) -c $(DBCSRMODS) $(DBCSR_MOD)

else

$(XMYMODULES): 		$(MYMODULES)
			$(FC) $(FFLAGS) $(PPFLAGS) -c $(MYMODULES)

endif

$(XMYSUBS): 		$(MYSUBS)
			$(FC) $(FFLAGS) $(PPFLAGS) -c $(MYSUBS)

clean:
			rm *.o *.mod *.a
