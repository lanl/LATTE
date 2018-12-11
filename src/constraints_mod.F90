
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


!> ! To add constrains to the system.
!! This module can be used to add constrains to the system.
!! For example, fixing the position of certain atoms.
!! \todo Add constrains for orientation and mass center distances.
!!
MODULE CONSTRAINTS_MOD

  USE CONSTANTS_MOD

  USE KERNELPARSER_MOD

#ifdef PROGRESSON
  USE bml
#endif

  IMPLICIT NONE

  PRIVATE

  INTEGER, PARAMETER :: DP = LATTEPREC

  PUBLIC :: FREEZE_ATOMS

  !> General input variables for constrains.
  !! \todo construct a parser for the constrains
  !!
  TYPE, PUBLIC :: CONSTRAINS_TYPE

     !> Verbosity level.
     INTEGER :: VERBOSE

     !> Constrain method
     CHARACTER(20) :: METHOD

  END TYPE CONSTRAINS_TYPE

CONTAINS

  !> To freeze a group of atoms.
  !! The atoms indices to be frozen are read from a file.
  !! The file is formated as follows:
  !!
  !! 3 #Number of freezed atoms
  !! 10  #We are freezing all the coordinates of atom 10
  !! 20  #We are freezing all the coordinates of atom 10
  !! 40  #We are freezing all the coordinates of atom 10
  !! \todo Generalize this to have constrains to different coordinates.
  !!
  !! \param FTOT Total forces. FTOT(1,3) gives the force on x direction for atom 3.
  !! \param VEL Velocities. VEL(1,3) gives the velocity on x direction for atom 3.
  !!
  SUBROUTINE FREEZE_ATOMS(FTOT,VEL)

    INTEGER, SAVE :: NFREEZE
    INTEGER, SAVE, ALLOCATABLE :: FREEZEID(:)
    INTEGER :: I
    REAL(DP), OPTIONAL, INTENT(INOUT) :: VEL(:,:)
    REAL(DP), INTENT(INOUT) :: FTOT(:,:)

    IF (.NOT. ALLOCATED(FREEZEID)) THEN
       OPEN(444,FILE="freeze.in")
       READ(444,*)NFREEZE
       ALLOCATE(FREEZEID(NFREEZE))
       DO I = 1,NFREEZE
          READ(444,*)FREEZEID(I)
       ENDDO
       CLOSE(444)
    ENDIF

    IF(PRESENT(VEL))THEN
       DO I = 1,NFREEZE
          VEL(:,FREEZEID(I)) = 0.0d0
          FTOT(:,FREEZEID(I)) = 0.0d0
       ENDDO
    ELSE
       DO I = 1,NFREEZE
          FTOT(:,FREEZEID(I)) = 0.0d0
       ENDDO
    ENDIF

  END SUBROUTINE FREEZE_ATOMS

END MODULE CONSTRAINTS_MOD
