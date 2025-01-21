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
SUBROUTINE LATTE_C_BIND(FLAGS, NATS, COORDS, TYPES, NTYPES, MASSES, XLO &
                        , XHI, XY, XZ, YZ, FORCES, MAXITER, VENERG, &
                        VEL, DT, VIRIAL_INOUT, NEWSYSTEM, EXISTERROR) BIND(C, NAME="latte")

  USE ISO_C_BINDING, ONLY: C_CHAR, C_NULL_CHAR, C_DOUBLE, C_INT, C_BOOL
  USE LATTE_LIB

  IMPLICIT NONE

  INTEGER(C_INT)                 ::  NATS, NTYPES, MAXITER
  INTEGER(C_INT)                 ::  TYPES(NATS), FLAGS(5)
  REAL(C_DOUBLE)                 ::  COORDS(3, NATS), MASSES(NTYPES), XHI(3)
  REAL(C_DOUBLE)                 ::  XLO(3), EKIN, VENERG, DT
  REAL(C_DOUBLE)                 ::  XY, XZ, YZ
  REAL(C_DOUBLE), INTENT(INOUT)  ::  FORCES(3, NATS), VEL(3, NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  ::  VIRIAL_INOUT(6)
  LOGICAL(C_BOOL)                ::  EXISTERROR
  INTEGER(C_INT), INTENT(INOUT)  ::  NEWSYSTEM

  CALL LATTE(NTYPES, TYPES, COORDS, MASSES, XLO, XHI, XY, XZ, YZ, FORCES, &
             MAXITER, VENERG, VEL, DT, VIRIAL_INOUT, NEWSYSTEM, EXISTERROR)

  RETURN

END SUBROUTINE LATTE_C_BIND

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
!! \brief This file is used to interface to python via iso_c_binding
!! library.
!! \param nats Number of total atoms in the system
!! \param nTypes Number of atom types
!! \param coords_in Coordinates of every atom in the system.
!! Allocation:
!! \verbatim coordinate(3,nats) \endverbatim
!! \param latticeVectors_in Flattened lattice vectors/box
!! Allocation:
!! \verbatim  lattice_vector(3*3) \endverbatim
!! \verbatim  v1 = lattice_vector(1:3) \endverbatim
!! \param atomTypes_in Atom type index for every atom
!! It gives the species index of a particulat atom. Indexing starts from 0!
!! Allocation:
!! \verbatim  atomTypes(nats) \endverbatim
!! If we need the index of atom 30 then:
!! \verbatim  atomTypes(30) \endverbatim
!! \param atomicNumbers_in Atomic number for every species.
!! A list with the atomic numbers for every species.
!! Allocation:
!! \verbatim  atomicNumbers(nTypes) \endverbatim
!! \return forces_out Computed forces. Flattened 2D array.
!! \verbatim forces_out(2) \endverbatim : y component of the force for atom 1.
!! \param charges_out Vector of computed charges.
!! \param verb_in Verbosity level.
!!
FUNCTION LATTE_COMPUTE(NATS, NTYPES, COORDS_IN, LATTICEVECTORS_IN, ATOMTYPES_IN, ATOMICNUMBERS_IN,&
     &FIELD_IN, CHARGES_OUT, FORCES_OUT, DIPOLE_OUT, ENERGY_OUT, VERB_IN) RESULT(ERR) BIND(C, NAME='latte_compute')
  USE LATTE_LIB
  USE ISO_C_BINDING, ONLY: C_CHAR, C_DOUBLE, C_INT, C_BOOL

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER(C_INT), INTENT(IN), VALUE  :: NATS
  INTEGER(C_INT), INTENT(IN), VALUE  :: NTYPES
  REAL(C_DOUBLE), INTENT(INOUT)  :: COORDS_IN(3*NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  :: FIELD_IN(3)
  REAL(C_DOUBLE), INTENT(INOUT)  :: FORCES_OUT(3*NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  :: CHARGES_OUT(NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  :: DIPOLE_OUT(3), ENERGY_OUT(1)
  INTEGER(C_INT), INTENT(INOUT)  :: ATOMTYPES_IN(NATS)
  INTEGER(C_INT), INTENT(INOUT) :: ATOMICNUMBERS_IN(NTYPES)
  REAL(C_DOUBLE), INTENT(INOUT) :: LATTICEVECTORS_IN(9)
  INTEGER(C_INT), INTENT(IN), VALUE :: VERB_IN
  LOGICAL(C_BOOL) :: ERR

  REAL(DP), ALLOCATABLE :: COORDS(:, :)
  REAL(DP), ALLOCATABLE :: FORCES(:, :), CHARGES(:), DIPOLE(:)
  !real(dp), allocatable :: bornch(:,:)
  REAL(DP), ALLOCATABLE :: FIELD(:)
  REAL(DP), ALLOCATABLE :: LATTICEVECTORS(:, :)
  INTEGER, ALLOCATABLE :: ATOMTYPES(:), ATOMICNUMBERS(:)
  INTEGER :: K
  INTEGER :: VERB
  LOGICAL :: ERR_STATUS

  ERR = .TRUE.
  ALLOCATE (COORDS(3, NATS))
  ALLOCATE (ATOMTYPES(NATS))
  ALLOCATE (ATOMICNUMBERS(NTYPES))
  ALLOCATE (LATTICEVECTORS(3, 3))
  ALLOCATE (CHARGES(NATS))
  ALLOCATE (FORCES(3, NATS))
  ALLOCATE (DIPOLE(3))
  ALLOCATE (FIELD(3))

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

  DO K = 1, NATS !We correct for the indexing
    ATOMTYPES(K) = ATOMTYPES_IN(K) + 1
  END DO

  FIELD = FIELD_IN

  VERB = VERB_IN

  CALL COMPUTE(COORDS, ATOMTYPES, ATOMICNUMBERS, LATTICEVECTORS,&
       &FIELD, CHARGES, FORCES, DIPOLE, ENERGY_OUT, VERB)

  !We vectorize/flatten the forces to send back to python
  DO K = 1, NATS
    FORCES_OUT((K - 1)*3 + 1) = FORCES(1, K)
    FORCES_OUT((K - 1)*3 + 2) = FORCES(2, K)
    FORCES_OUT((K - 1)*3 + 3) = FORCES(3, K)
  END DO

  CHARGES_OUT(:) = CHARGES(:)
  DIPOLE_OUT(:) = DIPOLE(:)

  DEALLOCATE (COORDS)
  DEALLOCATE (FORCES)
  DEALLOCATE (CHARGES)
  DEALLOCATE (LATTICEVECTORS)
  DEALLOCATE (ATOMTYPES)
  DEALLOCATE (ATOMICNUMBERS)
  DEALLOCATE (DIPOLE)

  ERR = ERR_STATUS

  RETURN

END FUNCTION LATTE_COMPUTE

!> Call the latte library from the python interface.
!! \brief This file is used to interface to python via iso_c_binding
!! library.
!! \param nats Number of total atoms in the system
!! \param nTypes Number of atom types
!! \param coords_in Coordinates of every atom in the system.
!! Allocation:
!! \verbatim coordinate(3,nats) \endverbatim
!! \param latticeVectors_in Flattened lattice vectors/box
!! Allocation:
!! \verbatim  lattice_vector(3*3) \endverbatim
!! \verbatim  v1 = lattice_vector(1:3) \endverbatim
!! \param atomTypes_in Atom type index for every atom
!! It gives the species index of a particulat atom. Indexing starts from 0!
!! Allocation:
!! \verbatim  atomTypes(nats) \endverbatim
!! If we need the index of atom 30 then:
!! \verbatim  atomTypes(30) \endverbatim
!! \param atomicNumbers_in Atomic number for every species.
!! A list with the atomic numbers for every species.
!! Allocation:
!! \verbatim  atomicNumbers(nTypes) \endverbatim
!! \return forces_out Computed forces. Flattened 2D array.
!! \verbatim forces_out(2) \endverbatim : y component of the force for atom 1.
!! \param charges_out Vector of computed charges.
!! \param verb_in Verbosity level.
!!
FUNCTION LATTE_COMPUTE_HS(NORBS, NATS, NTYPES, COORDS_IN, LATTICEVECTORS_IN, ATOMTYPES_IN, ATOMICNUMBERS_IN,&
     &FIELD_IN, HAM_OUT, OVER_OUT, VERB_IN) RESULT(ERR) BIND(C, NAME='latte_compute_hs')
  USE LATTE_LIB
  USE ISO_C_BINDING, ONLY: C_CHAR, C_DOUBLE, C_INT, C_BOOL

  IMPLICIT NONE
  INTEGER, PARAMETER :: DP = KIND(1.0D0)
  INTEGER(C_INT), INTENT(IN), VALUE  :: NORBS
  INTEGER(C_INT), INTENT(IN), VALUE  :: NATS
  INTEGER(C_INT), INTENT(IN), VALUE  :: NTYPES
  REAL(C_DOUBLE), INTENT(INOUT)  :: COORDS_IN(3*NATS)
  REAL(C_DOUBLE), INTENT(INOUT)  :: FIELD_IN(3)
  REAL(C_DOUBLE), INTENT(INOUT)  :: HAM_OUT(NORBS*NORBS)
  REAL(C_DOUBLE), INTENT(INOUT)  :: OVER_OUT(NORBS*NORBS)
  INTEGER(C_INT), INTENT(INOUT)  :: ATOMTYPES_IN(NATS)
  INTEGER(C_INT), INTENT(INOUT) :: ATOMICNUMBERS_IN(NTYPES)
  REAL(C_DOUBLE), INTENT(INOUT) :: LATTICEVECTORS_IN(9)
  INTEGER(C_INT), INTENT(IN), VALUE :: VERB_IN
  LOGICAL(C_BOOL) :: ERR

  REAL(DP), ALLOCATABLE :: COORDS(:, :)
  REAL(DP), ALLOCATABLE :: HAM(:, :), OVER(:, :)
  REAL(DP), ALLOCATABLE :: FIELD(:)
  REAL(DP), ALLOCATABLE :: LATTICEVECTORS(:, :)
  INTEGER, ALLOCATABLE :: ATOMTYPES(:), ATOMICNUMBERS(:)
  INTEGER :: K
  INTEGER :: VERB
  LOGICAL :: ERR_STATUS

  ERR = .TRUE.
  ALLOCATE (COORDS(3, NATS))
  ALLOCATE (ATOMTYPES(NATS))
  ALLOCATE (ATOMICNUMBERS(NTYPES))
  ALLOCATE (LATTICEVECTORS(3, 3))
  ALLOCATE (HAM(NORBS, NORBS))
  ALLOCATE (OVER(NORBS, NORBS))
  ALLOCATE (FIELD(3))

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

  DO K = 1, NATS !We correct for the indexing
    ATOMTYPES(K) = ATOMTYPES_IN(K) + 1
  END DO

  FIELD = FIELD_IN

  VERB = VERB_IN

  CALL COMPUTE_HS(NORBS, COORDS, ATOMTYPES, ATOMICNUMBERS, LATTICEVECTORS,&
       &FIELD, HAM, OVER, VERB)

  !We vectorize/flatten the forces to send back to python
  DO K = 1, NORBS
    HAM_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = HAM(:, K)
    OVER_OUT(1 + (K - 1)*NORBS:NORBS + (K - 1)*NORBS) = OVER(:, K)
  END DO
  !ham_out = reshape(ham, [size(ham)])
  !over_out = reshape(over, [size(over)])

  DEALLOCATE (COORDS)
  DEALLOCATE (HAM)
  DEALLOCATE (OVER)
  DEALLOCATE (LATTICEVECTORS)
  DEALLOCATE (ATOMTYPES)
  DEALLOCATE (ATOMICNUMBERS)
  DEALLOCATE (FIELD)

  ERR = ERR_STATUS

  RETURN

END FUNCTION LATTE_COMPUTE_HS
