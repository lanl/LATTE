#
# This is the main makefile for the serial version of LATTE
#
.SUFFIXES: .fortran90 .F

include ../makefile.CHOICES

PPFLAGS = -WF,-DXSYEVD,-DXHEEVD,-D$(PRECISION)PREC,-DGPU$(GPUOPT),-DDBCSR_$(DBCSR_OPT),-DMPI_$(MPIOPT)
PPFLAGS += -DPROGRESS$(PROGRESS)

# Flags for linking GPU libraries
GPU_LINKFLAG = -L../MATRIX -lxl 
GPU_LIB =  $(GPU_CUDA_LIB) -L../MATRIX -lmatrix_cuda.$(PRECISION) -lxlopt -lxlf90_r -lxlfmath -lxl -lxlsmp


MYMODULES = myprecision.F constants_mod.F setuparray.F coulombarray.F \
diagarray.F fermicommon.F mdarray.F neblistarray.F \
nonoarray.F ppotarray.F purearray.F relaxcommon.F restartarray.F \
sparsearray.F spinarray.F univarray.F virialarray.F xboarray.F \
kspacearray.F timer_mod.F subgraph.F matrixio.F sparsemath.F \
sparsesp2.F subgraphsp2.F homolumo.F openfiles_mod.F kernelparser_mod.F \
latteparser_latte_mod.F mixer_mod.F sp2progress.F genXprogress.F  \
bodirectprogress.F constraints_mod.F latte_lib.F 

DBCSRMODS = dbcsr_var_mod.F init_dbcsr.F shutdown_dbcsr.F

MYSUBS = allocatediag.F allocatenono.F allocatepure.F bodirect.F \
deallocatediag.F deallocatepure.F deorthomyrho.F deallocatenono.F \
diagmyh.F entropy.F fermiexpans.F genX.F  latte.F orthomyH.F \
 orthomyrho.F  pulay.F pulay_sp.F solvematcg.F solvematlapack.F \
sp2T.F pulay_spprogress.F \
sp2fermi.F sp2fermi_init.F sp2pure.F sp2pure_sparse.F kdiagmyh.F \
 kbodirect.F sp2gap.F sp2gap_setup.F readrestart.F initrng.F \
parawrite.F parafileopen.F wrtcfgs.F wrtrestart.F pbc.F \
addqdep.F allocatecoulomb.F allocatenebarrays.F  allocatexbo.F  \
atomcharge.F avepress.F avetemp.F bldnewH.F bldspinH.F \
conjgradient.F coulomb_ewald.F coulomb_oldskool.F coulomb_rspace.F \
coultailcoef.F deallocateall.F deallocatecoulomb.F \
deallocatenebarrays.F  deallocatexbo.F fcoulnono.F fcoulnono_sp.F \
fermiallocate.F fermideallocate.F \
fspinnono.F fspinnono_sp.F  gershgorin.F getbndfil.F \
getcoule.F getdensity.F getdeltaq.F getdeltaspin.F gethdim.F \
getke.F getmaxf.F getmdf.F getpressure.F getrho.F getspinE.F  \
 ifrestart.F initcoulomb.F initialv.F initshockcomp.F msrelax.F \
neblists.F noelec.F norms.F nptrescale.F  nvtrescale.F \
pairpot.F panic.F plot_ppot.F plot_univ.F printsparse.F \
 propchempot_xbo.F propspins_xbo.F qconsistency.F readcontrols.F \
readcr.F readmdcontroller.F readppot.F  readtb.F resetprodhd.F \
rhozero.F shockcomp.F  spinrhodirect.F stdescent.F \
summary.F tbmd.F toteng.F univscaling_function.F univtailcoef.F  \
vdwtailcoef.F velverlet.F  xbo.F am.F bm.F factorial.F \
 slmmp.F tlmmp.F wignerd.F dfda.F dfdb.F dfdr.F \
dunivscaling_function.F dslmmpda.F dslmmpdb.F dtlmmpda.F \
dtlmmpdb.F dwignerddb.F gradH.F kbldnewh.F  kgetrho.F kgradH.F  \
bldnewHS_sp.F gradH_sp.F univscaling.F dunivscaling.F kgetdos.F \
qneutral.F shiftH.F getmatindlist.F nvtlangevin.F nvtandersen.F \
nvtNH.F sp2pure_sparse_parallel.F sp2pure_sparse_parallel_simple.F \
sp2pure_subgraph_parallel.F gaussrn.F nnz.F allocatesubgraph.F \
deallocatesubgraph.F dosfit.F mofit.F getrespf.F gendiag.F \
mofit_plato.F ppfit.F pairpot_noneb.F allfit.F \
readppottab.F pairpottab.F hugrescale.F avesforhug.F getdipole.F \
fittingoutput.F tabtest.F kgenX.F korthomyH.F kdeorthomyrho.F \
kfcoulnono.F kpulay.F getforce.F masses2symbols.F setuptbmd.F \
wrtrestartlib.F deorthomyrhoprogress.F orthomyHprogress.F readrestartlib.F \
readppotspline.F pairpotspline.F errors.F initcoulombklist.F \
flcnnono.F flcnnono_sp.F kflcnnono.F dosfittab.F \
readppotD.F pairpot_plusD.F

# dependencies
# DEP_FILE = makefile.dep
# include $(DEP_FILE)

XMYMODULES = 		$(MYMODULES:.F=.o)
XDBCSRMODS = 		$(DBCSRMODS:.F=.o)
XMYSUBS = 		$(MYSUBS:.F=.o)


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
		$(FC) $(FFLAGS) -c latte_c_bind.F
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
