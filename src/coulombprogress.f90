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

!> To compute Coulombic contributions using routines from the progress library
!!
MODULE COULOMBPROGRESS_MOD

#ifdef PROGRESSON

  USE MYPRECISION

  PRIVATE

  PUBLIC :: COULOMBPRG

  LOGICAL, PUBLIC  ::  COULINIT = .FALSE.

CONTAINS

  !> This routine computes the real and reciprocal space contribution of the Ewald summation.
  !! \param atomi Atom index where Ewald Real will be calculated.
  !! \param spindex Species index list.
  !! \param splist Element symbol for every species.
  !! \param coordinates Coordinates for every atom in the system.
  !! \param charges Charges for every atom in the system.
  !! \param hubbardu Hubbard parameter U for every species.
  !! \param lattice_vectors Lattice vectors for the system.
  !! \param coul_acc Coulomb accuracy.
  !! \param coul_forces_k Coulombic forces (k space contribution)
  !! \param coul_pot_k Coulombic potential (k space contribution)
  !!
  SUBROUTINE COULOMBPRG(SPACE)
    IMPLICIT NONE
    CHARACTER(len=1), INTENT(IN) :: SPACE

    IF ( SPACE .EQ. "r" ) THEN
      CALL GET_EWALD_REAL(spindex,splist,coordinates,charges,hubbardu&
      &,lattice_vectors,volr,coul_acc,coul_forces_r,coul_pot_r)
    ELSE IF ( SPACE .EQ. "k") THEN
      CALL GET_EWALD_RECIP(spindex,splist,coordinates,charges,hubbardu,&
      &lattice_vectors,recip_vectors,volr,coul_acc,coul_forces_k,coul_pot_k)
    ENDIF

  END SUBROUTINE COULOMBPRG

#endif

END MODULE COULOMBPROGRESS_MOD
