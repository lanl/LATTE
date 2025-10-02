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

!> Subroutine for Latte-C++ interfacing.
!! \param FLAGS Different control flags that can be passed to LATTE (not in use yet)
!! \param NATS Number of atoms
!! \param COORDS Coordinates. Example: y-coordinate of atom 1 = COORDS(2,1)
!! \param TYPES An index for all the different atoms in the system.
!! \param NTYPES Number of different elements in the system
!! \param MASSES Element masses for every different element of the system.
!! \param XLO Lowest dimensions of the box
!! \param XHI Highest dimensions of the box
!! \param XY Tilt factor.
!! \param XZ Tilt factor.
!! \param YZ Tilt factor. The lattice vectors are constructed as:
!! a = (xhi-xlo,0,0); b = (xy,yhi-ylo,0); c = (xz,yz,zhi-zlo).
!! \param FORCES Forces for every atom as output.
!! \param MAXITER Latte MAXITER keyword. If MAXITER = -1, only the Forces are computed.
!!        If MAXITER = 0, MAXITER is read from latte.in file.
!!        IF MAXITER > 0, MAXITER is passed trough the library call.
!! \param VENERG This is the potential Energy that is given back from latte to the hosting code.
!! \param VEL Velocities passed to latte.
!! \param DT integration step passed to latte.
!! \param DT integration step passed to latte.
!! \param VIRIAL_INOUT Components of the second virial coefficient
!! \param NEWSYSTEM Tells LATTE if a new system is passed.
!! \param EXISTERROR Returns an error flag (.true.) to the hosting code.
!!
!! \brief This routine will be used load call latte_lib from a C/C++ program:
!!
!! \brief Note: To get the mass of atom 3 we do:
!! \verbatim
!!      MASS(TYPES(3))
!! \endverbatim
!!
!! \brief Note: To get the lattice vectors as formated in LATTE we do:
!! \verbatim
!!      BOX(1,1) = XHI(1) - XLO(1); ...
!! \endverbatim
!!
!! \brief Note: All units are LATTE units by default.
!! See https://github.com/losalamos/LATTE/blob/master/Manual/LATTE_manual.pdf
!!
FUNCTION LATTE_C_BIND(FLAGS_IN, NATS, COORDS_IN, TYPES_IN, NTYPES, MASSES_IN, XLO_IN &
                        , XHI_IN, XY, XZ, YZ, FORCES, MAXITER, VENERG_OUT, &
                        VEL_IN, DT, VIRIAL_INOUT, CHARGES_OUT, NEWSYSTEM_IN) RESULT(ERR) BIND(C, NAME="latte_c_bind")

  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_BOOL
  USE LATTE_LIB

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER(C_INT), INTENT(IN), VALUE  ::  NATS, NTYPES, MAXITER
  INTEGER(C_INT), INTENT(INOUT)  ::  TYPES_IN(NATS), FLAGS_IN(5)
  REAL(C_DOUBLE), INTENT(IN), VALUE  ::  XY, XZ, YZ, DT
  REAL(C_DOUBLE), INTENT(INOUT)  ::  COORDS_IN(3 * NATS), MASSES_IN(NTYPES), XHI_IN(3)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  XLO_IN(3), VENERG_OUT(1)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  FORCES(3, NATS), VEL_IN(3, NATS), CHARGES_OUT(NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  VIRIAL_INOUT(6)
  LOGICAL(C_BOOL) :: ERR  
  INTEGER(C_INT), INTENT(IN), VALUE  ::  NEWSYSTEM_IN
  
  REAL(DP) :: VENERG
  LOGICAL(1) :: ERR_STATUS
  INTEGER :: K

  INTEGER, ALLOCATABLE :: TYPES(:), FLAGS(:)
  REAL(DP), ALLOCATABLE :: COORDS(:, :), MASSES(:), XHI(:)
  REAL(DP), ALLOCATABLE :: XLO(:)
  REAL(DP), ALLOCATABLE :: VEL(:, :)

  ERR = .TRUE.
  ALLOCATE (TYPES(NATS))
  ALLOCATE (FLAGS(5))
  ALLOCATE (COORDS(3, NATS))
  ALLOCATE (MASSES(NTYPES))
  ALLOCATE (XHI(3))
  ALLOCATE (XLO(3))
  ALLOCATE (VEL(3, NATS))
  
  TYPES = TYPES_IN
  FLAGS = FLAGS_IN
!  COORDS = COORDS_IN
  MASSES = MASSES_IN
  XHI = XHI_IN
  XLO = XLO_IN
  VEL = VEL_IN

DO K = 1, NATS
 COORDS(1, K) = COORDS_IN((K - 1)*3 + 1)
 COORDS(2, K) = COORDS_IN((K - 1)*3 + 2)
 COORDS(3, K) = COORDS_IN((K - 1)*3 + 3)
END DO
!  write(*,*) NTYPES, TYPES, COORDS, MASSES, XLO, XHI, XY, XZ, YZ, FORCES, &
!       MAXITER, VENERG, VEL, DT, VIRIAL_INOUT, CHARGES_OUT, NEWSYSTEM_IN, ERR_STATUS
!  stop
  CALL LATTE(NTYPES, TYPES, COORDS, MASSES, XLO, XHI, XY, XZ, YZ, FORCES, &
       MAXITER, VENERG, VEL, DT, VIRIAL_INOUT, CHARGES_OUT, NEWSYSTEM_IN, ERR_STATUS)  
  VENERG_OUT = VENERG

  DEALLOCATE (TYPES)
  DEALLOCATE (FLAGS)
  DEALLOCATE (COORDS)
  DEALLOCATE (MASSES)
  DEALLOCATE (XHI)
  DEALLOCATE (XLO)
  DEALLOCATE (VEL)

  ERR = ERR_STATUS
 ! write(*,*) "VIRIAL_INOUT", VIRIAL_INOUT
  RETURN

END FUNCTION LATTE_C_BIND

!> Function for Latte-C++ interfacing.
!! \return ABIVERSION integer representing the date of the last change
!!          to the C/C++ interface (e.g. 20180221)
!!
!! \brief This function will be used prior to calling the LATTE library
!!        to allow the calling code to ensure the linked library version is compatible.
!!
INTEGER(C_INT) FUNCTION LATTE_C_ABIVERSION() BIND(C, NAME="latte_abiversion")
  USE ISO_C_BINDING, ONLY: C_INT
  USE LATTE_LIB, ONLY: LATTE_ABIVERSION
  IMPLICIT NONE

  LATTE_C_ABIVERSION = LATTE_ABIVERSION
  RETURN

END FUNCTION LATTE_C_ABIVERSION

!> Call the latte library from the python interface.
!! \param COMPFLAG_IN Different control flags that can be passed to LATTE
!!        COMPFLAG_IN = 1: Compute Hamiltonian and Overlap Matrices
!!        COMPFLAG_IN = 2: Compute 1 & Density Matrix, Charges, Eigenvalues, and Dvals 
!!        COMPFLAG_IN = 3: Compute 1 + 2 & Energy and Forces
!! \param NORBS Number atomic orbitals 
!! \param NCORES Number of cores in the graph partition scheme 
!! \param NATS Number of atoms
!! \param NTYPES Number of atom types
!! \param CHEMPOT_IN Chemical potential for the whole system
!! \param COULOMBV_IN Coulomb potential for each atom (NATS,) 
!! \param COORDS_IN Coordinates. Example: y-coordinate of atom 1 = COORDS_IN(2,1)
!! \param LATTICE_VECTORS_IN System PBC slab edge vectors. The first vector's coordinates 
!!        are LATTICE_VECTORS_IN(1,1:3) 
!! \param ATOMTYPES_IN An index for all the different atoms in the system (NATS,)
!! \param ATOMIC_NUMBERS_IN Atomic numbers for each species (NTYPES,)
!! \param HAM_OUT Hailtonian matrix as output (NORBS, NORBS)
!! \param OVER_OUT Overlap matrix as output (NORBS, NORBS)
!! \param DM_OUT Denisty matrix as output (NORBS, NORBS) 
!! \param CHARGES_OUT charges vector as output (NATS,)
!! \param EVALS_OUT Eigenvalues as output (NORBS,) 
!! \param DVALS_OUT Contribution to the eigenvectors of the system from the "core" as output (NORBS,) 
!! \param ENERGY_OUT This is the potential Energy that is given back from latte to the hosting code 
!! \param FORCES_OUT Forces for every atom as output (3, NATS)
!! \param VERB_IN Verbosity setting for LATTE 
!! \param NEWSYSTEM_IN Tells LATTE if a new system is passed.
!! \param EXISTERROR_INOUT Returns an error flag (.true.) to the hosting code.
!! \brief Note: All units are LATTE units by default. See https://github.com/losalamos/LATTE/blob/master/Manual/LATTE_manual.pdf
!! \brief Note: Control needs to be set to 1 to compute hamiltonian diagonalization when calling LATTE
!!
FUNCTION LATTE_COMPUTE(COMPFLAG_IN, SY_IDX, NUM_SY, NORBS, NCORES, NCOREATOMS, NATS, NTYPES, ETEMP, CHEMPOT_IN, COULOMBV_IN, COORDS_IN, LATTICEVECTORS_IN, ATOMTYPES_IN, ATOMICNUMBERS_IN,&
  & HAM_OUT, OVER_OUT, ZMAT_OUT, EVECTS_OUT, DM_OUT, CHARGES_OUT, EVALS_OUT, DVALS_OUT, ENERGY_OUT, FORCES_OUT, VERB_IN, NEWSYSTEM_IN, KEEPMEM_IN) RESULT(ERR) BIND(C, NAME='latte_compute')
USE LATTE_LIB
USE ISO_C_BINDING, ONLY: C_CHAR, C_DOUBLE, C_INT, C_BOOL
! USE NVTX_MOD

IMPLICIT NONE
INTEGER, PARAMETER :: DP = KIND(1.0D0)
INTEGER(C_INT), INTENT(IN), VALUE  :: COMPFLAG_IN
INTEGER(C_INT), INTENT(IN), VALUE  :: NATS
INTEGER(C_INT), INTENT(IN), VALUE  :: SY_IDX, NUM_SY 
INTEGER(C_INT), INTENT(IN), VALUE  :: NORBS
INTEGER(C_INT), INTENT(IN), VALUE  :: NTYPES
REAL(C_DOUBLE), INTENT(IN), VALUE  :: ETEMP 
REAL(C_DOUBLE), INTENT(IN), VALUE  :: CHEMPOT_IN 
REAL(C_DOUBLE), INTENT(INOUT)  :: COULOMBV_IN(NATS)
REAL(C_DOUBLE), INTENT(INOUT)  :: COORDS_IN(3*NATS)
REAL(C_DOUBLE), INTENT(INOUT)  :: FORCES_OUT(3*NATS)
REAL(C_DOUBLE), INTENT(INOUT)  :: CHARGES_OUT(NATS)
REAL(C_DOUBLE), INTENT(INOUT)  :: ENERGY_OUT(1)
REAL(C_DOUBLE), INTENT(INOUT)  :: HAM_OUT(NORBS*NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: OVER_OUT(NORBS*NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: ZMAT_OUT(NORBS*NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: EVECTS_OUT(NORBS*NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: DM_OUT(NORBS*NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: EVALS_OUT(NORBS)
REAL(C_DOUBLE), INTENT(INOUT)  :: DVALS_OUT(NORBS)
INTEGER(C_INT), INTENT(INOUT)  :: ATOMTYPES_IN(NATS)
INTEGER(C_INT), INTENT(INOUT) :: ATOMICNUMBERS_IN(NTYPES)
REAL(C_DOUBLE), INTENT(INOUT) :: LATTICEVECTORS_IN(9)
INTEGER(C_INT), INTENT(IN), VALUE :: VERB_IN, NEWSYSTEM_IN, KEEPMEM_IN
INTEGER(C_INT), INTENT(IN), VALUE :: NCORES, NCOREATOMS
LOGICAL(C_BOOL) :: ERR

REAL(DP), ALLOCATABLE :: COULOMBV(:)
REAL(DP), ALLOCATABLE :: COORDS(:, :)
REAL(DP), ALLOCATABLE :: FORCES(:, :), CHARGES(:)
REAL(DP), ALLOCATABLE :: HAM(:, :), OVER(:,:), ZMAT(:,:), EVECTS(:,:), DM(:,:)
REAL(DP), ALLOCATABLE :: EVALS(:), DVALS(:)
REAL(DP), ALLOCATABLE :: LATTICEVECTORS(:, :)
INTEGER, ALLOCATABLE :: ATOMTYPES(:), ATOMICNUMBERS(:)
INTEGER :: K
REAL(DP), ALLOCATABLE :: ENERGY(:)
LOGICAL(1) :: ERR_STATUS

ERR = .TRUE.
ALLOCATE (COULOMBV(NATS))
ALLOCATE (COORDS(3, NATS))
ALLOCATE (ATOMTYPES(NATS))
ALLOCATE (ATOMICNUMBERS(NTYPES))
ALLOCATE (LATTICEVECTORS(3, 3))
ALLOCATE (CHARGES(NATS))
ALLOCATE (ENERGY(1))
ALLOCATE (FORCES(3, NATS))
ALLOCATE (HAM(NORBS, NORBS))
ALLOCATE (OVER(NORBS, NORBS))
ALLOCATE (ZMAT(NORBS, NORBS))
ALLOCATE (EVECTS(NORBS, NORBS))
ALLOCATE (DM(NORBS, NORBS))
ALLOCATE (EVALS(NORBS))
ALLOCATE (DVALS(NORBS))

!Note that arrays appear in another order. We need to rearange
!the data. This is because of the column mayor (in python) vs.
!row mayor in fortran.
DO K = 1, NATS
 COORDS(1, K) = COORDS_IN((K - 1)*3 + 1)
 COORDS(2, K) = COORDS_IN((K - 1)*3 + 2)
 COORDS(3, K) = COORDS_IN((K - 1)*3 + 3)
END DO

LATTICEVECTORS(1, 1) = LATTICEVECTORS_IN(1)
LATTICEVECTORS(1, 2) = LATTICEVECTORS_IN(2)
LATTICEVECTORS(1, 3) = LATTICEVECTORS_IN(3)

LATTICEVECTORS(2, 1) = LATTICEVECTORS_IN(4)
LATTICEVECTORS(2, 2) = LATTICEVECTORS_IN(5)
LATTICEVECTORS(2, 3) = LATTICEVECTORS_IN(6)

LATTICEVECTORS(3, 1) = LATTICEVECTORS_IN(7)
LATTICEVECTORS(3, 2) = LATTICEVECTORS_IN(8)
LATTICEVECTORS(3, 3) = LATTICEVECTORS_IN(9)

ATOMICNUMBERS = ATOMICNUMBERS_IN
COULOMBV = COULOMBV_IN

DO K = 1, NATS !We correct for the indexing
 ATOMTYPES(K) = ATOMTYPES_IN(K) + 1
END DO
! call nvtxStartRange("LATTE_COMPUTE",1)
CALL COMPUTE(COMPFLAG_IN, SY_IDX, NUM_SY, NCORES, NCOREATOMS, ETEMP, CHEMPOT_IN, COULOMBV, COORDS, ATOMTYPES, ATOMICNUMBERS, LATTICEVECTORS,&
    & HAM, OVER, ZMAT, EVECTS, DM, CHARGES, EVALS, DVALS, ENERGY, FORCES, VERB_IN, NEWSYSTEM_IN, KEEPMEM_IN, ERR_STATUS)
! call nvtxEndRange
!We vectorize/flatten the forces to send back to python
DO K = 1, NATS
 FORCES_OUT((K - 1)*3 + 1) = FORCES(1, K)
 FORCES_OUT((K - 1)*3 + 2) = FORCES(2, K)
 FORCES_OUT((K - 1)*3 + 3) = FORCES(3, K)
END DO
DO K = 1, NORBS
  HAM_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = HAM(:, K)
  OVER_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = OVER(:, K)
  ZMAT_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = ZMAT(:, K)
  EVECTS_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = EVECTS(:, K)
  DM_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = DM(:, K)
END DO

CHARGES_OUT(:) = CHARGES(:)
EVALS_OUT(:) = EVALS(:)
DVALS_OUT(:) = DVALS(:)
ENERGY_OUT(:) = ENERGY(:)

DEALLOCATE (COULOMBV)
DEALLOCATE (COORDS)
DEALLOCATE (LATTICEVECTORS)
DEALLOCATE (ATOMTYPES)
DEALLOCATE (ATOMICNUMBERS)
DEALLOCATE (CHARGES)
DEALLOCATE (ENERGY)
DEALLOCATE (FORCES)
DEALLOCATE (HAM)
DEALLOCATE (OVER)
DEALLOCATE (ZMAT)
DEALLOCATE (EVECTS)
DEALLOCATE (DM)
DEALLOCATE (EVALS)
DEALLOCATE (DVALS)

ERR = ERR_STATUS

RETURN

END FUNCTION LATTE_COMPUTE
