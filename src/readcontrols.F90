  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Copyright 2010.  Los Alamos National Security, LLC. This material was    !
! produced under U.S. Government contract DE-AC52-06NA25396 for Los Alamos !
! National Laboratory (LANL), which is operated by Los Alamos National     !
! Security, LLC for the U.S. Department of Energy. The U.S. Government has !
! rights to use, reproduce, and distribute this software.  NEITHER THE     !
! GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY WARRANTY,     !
! EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS         !
! SOFTWARE.  If software is modified to produce derivative works, such     !
! modified software should be clearly marked, so as not to confuse it      !
! with the version available from LANL.                                    !
!                                                                          !
! Additionally, this program is free software; you can redistribute it     !
! and/or modify it under the terms of the GNU General Public License as    !
! published by the Free Software Foundation; version 2.0 of the License.   !
! Accordingly, this program is distributed in the hope that it will be     !
! useful, but WITHOUT ANY WARRANTY; without even the implied warranty of   !
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General !
! Public License for more details.                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE READCONTROLS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE COULOMBARRAY
  USE FERMICOMMON
  USE SPARSEARRAY
  USE RELAXCOMMON

  IMPLICIT NONE

  CHARACTER(LEN=20) :: HD

  IF (EXISTERROR) RETURN

  OPEN(UNIT=13, STATUS="OLD", FILE=TRIM(PARAMPATH)//"/control.in")

  !
  ! CONTROL determines how the density matrix is going to be
  ! calculated: 1 = diagonalization, 2 = SP2 purification,
  ! 3 = recursive expansion of the Fermi operator, 4 = SP2T,
  ! 5 = SP2 Fermi (truncated SP2)
  !

  READ(13,*) HD, CONTROL

  !
  ! BASISTYPE can equal "ORTHO" OR "NONORTHO",
  !

  READ(13,*) HD, BASISTYPE

  IF (BASISTYPE .NE. "ORTHO" .AND. BASISTYPE .NE. "NONORTHO") THEN
     CALL ERRORS("readcontrols","Error defining basis type (ortho/nonortho)")
  ENDIF

  READ(13,*) HD, SCLTYPE ! Choose whether to use the analytic forms or tables

  READ(13,*) HD, DEBUGON


  !
  ! Read the order of the recursion in the expansion of the Fermi
  ! operator, M.
  !

  READ(13,*) HD, FERMIM

  ! If we're using the expansion of the Fermi operator, we can
  ! use a LAPACK routine or Niklasson's conjugate gradient method to
  ! solve AX = B. CGORLIB: 0 = LAPACK, 1 = conjugate gradient
  ! CGTOL = the user-supplied tolerance for the CG solution of AX = B

  READ(13,*) HD, CGORLIB, HD, CGTOL

  CGTOL2 = CGTOL*CGTOL

  ! Electronic temperature, in eV

  READ(13,*) HD, KBT

  !
  ! Read the number of recursions for the truncated, finite
  ! temperature SP2 algorithm
  !

  READ(13,*) HD, NORECS

  !
  ! What kind of entropy are we going to use in a finite Te calculation
  !
  ! ENTROPYKIND = 0 : none
  ! ENTROPYKIND = 1 : exact for Fermi-Dirac occupation
  ! ENTROPYKIND = 2 : Different form of exact expression that may be useful
  ! when using CONTROL = 5
  ! ENTROPYKIND = 3 : 4th order expansion of exact form (no diag)
  ! ENTROPYKIND = 4 : 8th order expansion of exact form (no diag)
  !

  READ(13,*) HD, ENTROPYKIND

  !
  ! Do we want long-range C/R^6 tails?
  !
  ! PPOTON = 1: Turn on pairwise interaction
  ! PPOTON = 0: Turn it off (useful for fitting)
  !
  ! VDWON = 0: No C/R^6 tails
  ! VDWON = 1: Use tails
  !

  READ(13,*) HD, PPOTON, HD, VDWON


  !
  ! Are we doing a spin-polarized calculation?
  ! SPINON = 1 = yes
  ! SPINON = 0 = no

  READ(13,*) HD, SPINON, HD, SPINTOL

  !
  ! Controls for electrostatics:
  !
  ! ELECTRO: 0 = LCN is applied, 1 = charge dependent TB on
  ! ELECMETH: 0 = Ewald  summation, 1 = All real space
  ! ELEC_ETOL: Tolerance on energy when determining charges (not implemented)
  ! ELEC_QTOL: Tolerance on charges during self-consistent calc
  !

  READ(13,*) HD, ELECTRO, HD, ELECMETH, HD, ELEC_ETOL, HD, ELEC_QTOL

  !
  ! COULACC: Accuracy for the Ewald method (1.0e-4 works)
  ! COULCUT: If we're using the Ewald method, this is the cut-off for the
  ! real space part. If we're doing it all in real space, this is the radial
  ! cut-off for the sum.
  ! COULR1: If we're doing it in real space, the cut-off tail on 1/R is
  ! applied here at terminated at COULCUT.
  !

  READ(13,*) HD, COULACC, HD, COULCUT, HD, COULR1

  !
  ! MAXSCF:  Maximum number of SCF cycles
  !

  READ(13,*) HD, MAXSCF

  !
  ! BREAKTOL: Tolerance for breaking SP2 loops
  ! MINSP2ITER: Minimum number of iterations during SP2 purification
  !

  READ(13,*) HD, BREAKTOL, HD, MINSP2ITER, HD, SP2CONV

  !
  ! FULLQCONV: 0 = We'll run QITER SCF cycles during MD, = 1, we'll run
  ! SCF cycles until we've reached ELEC_QTOL. Only important for MD
  ! QITER: Number of SCF cycles we're going to run at each MD time step
  !

  READ(13,*) HD, FULLQCONV, HD, QITER

  !
  ! QMIX AND SPINMIX are the coefficients for the linear mixing of
  ! new and old charge and spin densities, respectively, during SCF cycles
  !

  READ(13,*) HD, QMIX, HD, SPINMIX, HD, MDMIX

  !
  ! ORDERNMOL: Turn on molecule-ID-based density matrix blocking
  !

  READ(13,*) HD, ORDERNMOL

  !
  ! SPARSEON: 0 = all dense matrix stuff, 1 = use CSR format and
  ! Gustavson's algorithm for matrix-matrix multiplication
  ! THRESHOLDON: 0 = do not throw away elements; 1 = throw away elements
  ! NUMTHRESH: If THRESHOLDON = 1 throw away element whose absolute value is
  ! smaller than NUMTHRESH
  ! FILLINSTOP: Number of purification cycles beyond which we stop allowing
  ! for further fill-in
  !

  READ(13,*) HD, SPARSEON, HD, THRESHOLDON,  HD, NUMTHRESH, HD, FILLINSTOP, HD, BLKSZ

  !
  ! MSPARSE: value for M when SPARSEON = 1, used by sp2 sparse algorithm
  !          0 = value for M is not known, defaults to N
  !

  READ(13,*) HD, MSPARSE

  !
  ! LCNON: 0 = during charge neutral MD simulations we'll run LCNITER SCF
  ! cycles at each time step, 1 = we'll run SCF cycles until CHTOL is reached
  ! LCNITER: Number of SCF cycles to achieve LCN at each MD time step
  ! CHTOL: Tolerance on atomic charges (Mulliken) before LCN is declared
  !

  READ(13,*) HD, LCNON, HD, LCNITER, HD, CHTOL

  !
  ! Read the SKIN for the neighbor list (Angstrom)
  !

  READ(13,*) HD, SKIN

  !
  ! RELAXME: 0 = Don't run relaxation, 1 = relax geometry
  ! RELTYPE: SD = steepest descent, CG = conjugate gradient
  ! MXRLX: Maximum number of steps in the geometry optimization
  ! RLXFTOT: Run optimization until all forces are less than RLXFTOL
  !

  READ(13,*) HD, RELAXME, HD, RELTYPE, HD, MXRLX, HD, RLXFTOL

  !
  ! MDON: 0 = Molecular dynamics off, 1 = Molecular dynamics on
  ! (MD is controlled using the file MDcontroller)
  !

  READ(13,*) HD, MDON

  !
  ! PBCON: 1 = full periodic boundary conditions, 0 = gas phase: no pbc and
  ! electrostatics done all in real space
  !

  READ(13,*) HD, PBCON

  READ(13,*) HD, RESTART

  ! Add or remove electrons. 2+ -> charge = +2 since TOTNE = TOTNE - CHARGE

  READ(13,*) HD, CHARGE

  !
  ! XBOON: 0 = Niklasson's extended Lagrangian Born-Oppenheimer MD off,
  ! 1 = on.
  !

  READ(13,*) HD, XBOON

  !
  ! XBODISON: We have the option of turning on damping for the XBO
  ! to remedy the accumulation of noise. 0 = off, 1 = on.
  !

  READ(13,*) HD, XBODISON

  !
  ! XBODISORDER: = Order of the damping function (1 - 9)
  !

  READ(13,*) HD, XBODISORDER

  FITON = 0

  !
  ! Read in the number of GPUs per node
  !

  READ(13,*) HD, NGPU

  ! Are we doing k-space?

  READ(13,*) HD, KON

  ! Do we want to calculate forces too (not always necessary when fitting)

  READ(13,*) HD, COMPFORCE

  ! Turn on the simulated annealing subroutine to fit DOS


  READ(13,*) HD, DOSFITON, HD, INT2FIT, HD, MCBETA, HD, NFITSTEP, HD, QFIT, &
       HD, MCSIGMA

  READ(13,*) HD, PPFITON

  READ(13,*) HD, ALLFITON

  READ(13,*) HD, PPNFITSTEP, HD, BINFITSTEP, HD, PP2FIT, HD, BINT2FIT

  READ(13,*) HD, PPBETA, HD, PPSIGMA, HD, PPNMOL, HD, PPNGEOM

  READ(13,*) HD, PARREP

  ! Dielectric constant

  READ(13,*) HD, RELPERM

  CLOSE(13)

  !
  ! Summarize the calculation we're doing here
  !

  !  OPEN(UNIT=99, STATUS="UNKNOWN", FILE="my_last_LATTE_calc")

  !  IF (CONTROL .EQ. 1) THEN
  !     WRITE(99,*) "Diagonization used to calculate bond order"
  !  ELSEIF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 0) THEN
  !     WRITE(99,*) "Dense matrix SP2 used to calculate bond order"
  !  ELSEIF (CONTROL .EQ. 2 .AND. SPARSEON .EQ. 1) THEN
  !     WRITE(99,*) "Quasi-sparse matrix SP2 used to calculated bond order"
  !  ELSEIF (CONTROL .EQ. 3) THEN
  !     WRITE(99,*) "Recursive expansion of the Fermi operator"
  !     IF (CGORLIB .EQ. 0) THEN
  !        WRITE(99,*) "Dense matrix: using LAPACK routine to solve AX=B"
  !     ENDIF
  !     IF (CGORLIB .EQ. 1 .AND. SPARSEON .EQ. 0) THEN
  !        WRITE(99,*) "Dense matrix conjugate gradient scheme to solve AX=B"
  !     ELSEIF (CGORLIB .EQ. 1 .AND. SPARSEON .EQ. 1) THEN
  !        WRITE(99,*) "Sparse matrix conjugate gradient scheme to solve AX=B"
  !     ENDIF
  !  ENDIF

  RETURN

END SUBROUTINE READCONTROLS
