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

MODULE CONSTANTS_MOD

  USE MYPRECISION

  IMPLICIT NONE
  SAVE

  INTEGER :: NATS, NMAT, NOINT, NOELEM, HDIM
  INTEGER :: CONTROL, RELAXME, PBCON, RESTART, MDON
  INTEGER :: NRANK  !! ANDERS_CHANGE
  INTEGER :: XBOON, XBODISON, XBODISORDER
  INTEGER :: SCFS, SCFS_II
  INTEGER :: LCNON, LCNITER
  INTEGER :: ELECTRO, ELECMETH
  INTEGER :: PPOTON, PLUSDON, VDWON
  INTEGER :: FULLQCONV, QITER
  INTEGER :: EMETH, CGORLIB
  INTEGER :: SPINON, SPARSEON
  INTEGER :: NSPIN
  INTEGER :: ORDERNMOL
  INTEGER :: MAXSCF, MINSP2ITER
  INTEGER :: VARDT
  INTEGER :: NORECS, ENTROPYKIND
  INTEGER :: DEBUGON, FITON
  INTEGER :: BLKSZ ! Block size when using DBCSR
  INTEGER :: NGPU ! Number of GPUs (blimey!)
  INTEGER :: NUMSCF
  INTEGER :: CHARGE
  INTEGER :: KON ! K-SPACE FLAG
  INTEGER :: COMPFORCE
  INTEGER :: SPONLY ! If we only have sp-bonded elements: faster gradH
  INTEGER :: DOSFITON, INT2FIT, NFITSTEP, QFIT ! For the simulated annealing
  INTEGER :: PP2FIT, BINT2FIT, PPNFITSTEP, PPNMOL, PPNGEOM, BINFITSTEP
  INTEGER :: PPFITON, ALLFITON
  INTEGER :: MDADAPT, MDADAPT_COUNT 
  INTEGER :: PARREP
  INTEGER :: VERBOSE = 0
  INTEGER :: MIXER
  INTEGER :: RSLEVEL
  INTEGER :: RESTARTLIB
  INTEGER :: FREEZE, FIRE_FREEZE
  INTEGER :: NEWSYSTEMLATTE 
  INTEGER :: KERNELSCHEME
  REAL(LATTEPREC) :: BOX(3,3), BOX_OLD(3,3), BOXDIMS(3)
  REAL(LATTEPREC) :: BNDFIL, TOTNE, MAXDN2DT
  REAL(LATTEPREC) :: COVE, TOTE, ENTE, KEE, ECOUL, EREP, TRRHOH
  REAL(LATTEPREC) :: TRRHOH0  ! ANDERS CHANGE
  REAL(LATTEPREC) :: ESPIN, ESPIN_ZERO
  REAL(LATTEPREC) :: EPLUSD, RESIDUEOLD
  REAL(LATTEPREC) :: MINEVAL, MAXEVAL, MAXMINUSMIN
  REAL(LATTEPREC) :: CHEMPOT, KBT
  REAL(LATTEPREC) :: EGAP, EHOMO, ELUMO, RESNORM, RESNORMOLD
  REAL(LATTEPREC) :: CHTOL, SPINTOL
  REAL(LATTEPREC) :: ELEC_ETOL, ELEC_QTOL
  REAL(LATTEPREC) :: QMIX, SPINMIX, MDMIX
  REAL(LATTEPREC) :: BREAKTOL
  REAL(LATTEPREC) :: SUMMASS, MASSDEN
  REAL(LATTEPREC) :: MCBETA, MCSIGMA ! Temperature in simulated annealing
  REAL(LATTEPREC) :: PPBETA, PPSIGMA
  REAL(LATTEPREC) :: MINGAP ! Minimum gap for adaptive SCF
  REAL(LATTEPREC) :: KERNELTOL ! tolerance for adaptive kernels
  REAL(LATTEPREC) :: S_DFTB_U ! S_DFTB_U is the chosen Hubbard U for the DM term
  REAL(LATTEPREC) :: P_DFTB_U ! P_DFTB_U is the chosen Hubbard U for the DM term
  REAL(LATTEPREC) :: D_DFTB_U ! D_DFTB_U is the chosen Hubbard U for the DM term
  REAL(LATTEPREC) :: F_DFTB_U ! F_DFTB_U is the chosen Hubbard U for the DM term
  CHARACTER(LEN = 3) :: SP2CONV
  CHARACTER(LEN = 20) :: BASISTYPE, SCLTYPE !scltype controls whether we do exp or tabular integrals

  CHARACTER(LEN = 100) :: PARAMPATH = "./TBparam"
  CHARACTER(LEN = 100) :: COORDSFILE = "./bl/inputblock.dat"

  ! For the latte lib
  CHARACTER(LEN = 20) :: JOB
  LOGICAL :: LIBINIT = .FALSE.
  INTEGER :: LIBCALLS = 0
  LOGICAL(1) :: EXISTERROR
  LOGICAL :: LIBRUN = .FALSE.
  LOGICAL :: STOPATMAXSCF

  !For the new input file parser
  CHARACTER(LEN = 1000) :: LATTEINNAME = "latte.in"
  CHARACTER(LEN = 1000) :: OUTFILE = "log.latte"
  LOGICAL :: LATTEINEXISTS
  LOGICAL :: DOKERNEL
  LOGICAL :: DFTBU
  LOGICAL :: READKERNEL
  LOGICAL :: SAVEKERNEL

  !For failsafe mode
  LOGICAL :: FAILSAFE

  !! For truncated SP2 and entropy calculation
  INTEGER :: SCFSTEP = 0
  INTEGER :: OCCSTEPS
  REAL(LATTEPREC) :: TSCALE
  REAL(LATTEPREC) :: OCCERRLIMIT, TRACELIMIT, EPS

  ! Some often-used constants
  REAL(LATTEPREC), PARAMETER :: MVV2KE = 166.0538782D0/1.602176487D0
  REAL(LATTEPREC), PARAMETER :: KE2T = 1.0D0/0.000086173435D0
  REAL(LATTEPREC), PARAMETER :: F2V = 9.6484504393669415d-003
  REAL(LATTEPREC), PARAMETER :: TOGPA = 160.2176487D0
  REAL(LATTEPREC), PARAMETER :: MVV2T = 1660538.782D0/1.3806504D0
  !  REAL(LATTEPREC), PARAMETER :: PI = 3.14159265358979323846264D0
  REAL(LATTEPREC), PARAMETER :: PI = TWO*ACOS(ZERO)


END MODULE CONSTANTS_MOD
