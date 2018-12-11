#
# This is the main makefile for the serial version of LATTE
#
.SUFFIXES: .fortran90 .F90

include ../makefile.CHOICES

PPFLAGS = -DXSYEVD -D$(PRECISION)PREC -DGPU$(GPUOPT) -DDBCSR_$(DBCSR_OPT) -DMPI_$(MPIOPT)
PPFLAGS += -DPROGRESS$(PROGRESS)

# Flags for linking GPU libraries
GPU_LINKFLAG = -L../MATRIX -lgfortran -lstdc++
GPU_LIB = -L../MATRIX -lmatrix_cuda.$(PRECISION) $(GPU_CUDA_LIB)


MYMODULES = myprecision.F90 constants_mod.F90 setuparray.F90 coulombarray.F90 \
diagarray.F90 fermicommon.F90 mdarray.F90 neblistarray.F90 \
nonoarray.F90 ppotarray.F90 purearray.F90 relaxcommon.F90 restartarray.F90 \
sparsearray.F90 spinarray.F90 univarray.F90 virialarray.F90 xboarray.F90 \
kspacearray.F90 timer_mod.F90 subgraph.F90 matrixio.F90 sparsemath.F90 \
sparsesp2.F90 subgraphsp2.F90 homolumo.F90 openfiles_mod.F90 kernelparser_mod.F90 \
latteparser_latte_mod.F90 mixer_mod.F90 sp2progress.F90 genXprogress.F90  \
bodirectprogress.F90 constraints_mod.F90 latte_lib.F90 

DBCSRMODS = dbcsr_var_mod.F90 init_dbcsr.F90 shutdown_dbcsr.F90

MYSUBS = allocatediag.F90 allocatenono.F90 allocatepure.F90 bodirect.F90 \
deallocatediag.F90 deallocatepure.F90 deorthomyrho.F90 deallocatenono.F90 \
diagmyh.F90 entropy.F90 fermiexpans.F90 genX.F90  latte.F90 orthomyH.F90 \
 orthomyrho.F90  pulay.F90 pulay_sp.F90 solvematcg.F90 solvematlapack.F90 \
sp2T.F90 pulay_spprogress.F90 \
sp2fermi.F90 sp2fermi_init.F90 sp2pure.F90 sp2pure_sparse.F90 kdiagmyh.F90 \
 kbodirect.F90 sp2gap.F90 sp2gap_setup.F90 readrestart.F90 initrng.F90 \
parawrite.F90 parafileopen.F90 wrtcfgs.F90 wrtrestart.F90 pbc.F90 \
addqdep.F90 allocatecoulomb.F90 allocatenebarrays.F90  allocatexbo.F90  \
atomcharge.F90 avepress.F90 avetemp.F90 bldnewH.F90 bldspinH.F90 \
conjgradient.F90 coulomb_ewald.F90 coulomb_oldskool.F90 coulomb_rspace.F90 \
coultailcoef.F90 deallocateall.F90 deallocatecoulomb.F90 \
deallocatenebarrays.F90  deallocatexbo.F90 fcoulnono.F90 fcoulnono_sp.F90 \
fermiallocate.F90 fermideallocate.F90 \
fspinnono.F90 fspinnono_sp.F90  gershgorin.F90 getbndfil.F90 \
getcoule.F90 getdensity.F90 getdeltaq.F90 getdeltaspin.F90 gethdim.F90 \
getke.F90 getmaxf.F90 getmdf.F90 getpressure.F90 getrho.F90 getspinE.F90  \
 ifrestart.F90 initcoulomb.F90 initialv.F90 initshockcomp.F90 msrelax.F90 \
neblists.F90 noelec.F90 norms.F90 nptrescale.F90  nvtrescale.F90 \
pairpot.F90 panic.F90 plot_ppot.F90 plot_univ.F90 printsparse.F90 \
 propchempot_xbo.F90 propspins_xbo.F90 qconsistency.F90 readcontrols.F90 \
readcr.F90 readmdcontroller.F90 readppot.F90  readtb.F90 resetprodhd.F90 \
rhozero.F90 shockcomp.F90  spinrhodirect.F90 stdescent.F90 \
summary.F90 tbmd.F90 toteng.F90 univscaling_function.F90 univtailcoef.F90  \
vdwtailcoef.F90 velverlet.F90  xbo.F90 am.F90 bm.F90 factorial.F90 \
 slmmp.F90 tlmmp.F90 wignerd.F90 dfda.F90 dfdb.F90 dfdr.F90 \
dunivscaling_function.F90 dslmmpda.F90 dslmmpdb.F90 dtlmmpda.F90 \
dtlmmpdb.F90 dwignerddb.F90 gradH.F90 kbldnewh.F90  kgetrho.F90 kgradH.F90  \
bldnewHS_sp.F90 gradH_sp.F90 univscaling.F90 dunivscaling.F90 kgetdos.F90 \
qneutral.F90 shiftH.F90 getmatindlist.F90 nvtlangevin.F90 nvtandersen.F90 \
nvtNH.F90 sp2pure_sparse_parallel.F90 sp2pure_sparse_parallel_simple.F90 \
sp2pure_subgraph_parallel.F90 gaussrn.F90 nnz.F90 allocatesubgraph.F90 \
deallocatesubgraph.F90 dosfit.F90 mofit.F90 getrespf.F90 gendiag.F90 \
mofit_plato.F90 ppfit.F90 pairpot_noneb.F90 allfit.F90 \
readppottab.F90 pairpottab.F90 hugrescale.F90 avesforhug.F90 getdipole.F90 \
fittingoutput.F90 tabtest.F90 kgenX.F90 korthomyH.F90 kdeorthomyrho.F90 \
kfcoulnono.F90 kpulay.F90 getforce.F90 masses2symbols.F90 setuptbmd.F90 \
wrtrestartlib.F90 deorthomyrhoprogress.F90 orthomyHprogress.F90 readrestartlib.F90 \
readppotspline.F90 pairpotspline.F90 errors.F90 initcoulombklist.F90 \
flcnnono.F90 flcnnono_sp.F90 kflcnnono.F90 dosfittab.F90

# dependencies
# DEP_FILE = makefile.dep
# include $(DEP_FILE)

XMYMODULES = 		$(MYMODULES:.F90=.o)
XDBCSRMODS = 		$(DBCSRMODS:.F90=.o)
XMYSUBS = 		$(MYSUBS:.F90=.o)


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
		$(FC) $(FFLAGS) -c latte_c_bind.F90
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
