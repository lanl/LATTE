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

MODULE FERMICOMMON

  USE MYPRECISION

  IMPLICIT NONE

  SAVE

  INTEGER :: FERMIM
  REAL(LATTEPREC) :: CGTOL, CGTOL2

  ! These are for the dense matrix version...

  !  REAL(LATTEPREC), ALLOCATABLE :: X(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: R0(:,:), P0(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: A(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: TMPMAT(:,:), X2(:,:)

  ! And these are for the sparse matrix one

  REAL(LATTEPREC), ALLOCATABLE :: VALRHO(:), VALR0(:), VALP0(:)
  REAL(LATTEPREC), ALLOCATABLE :: VALA(:), VALTMP(:)

END MODULE FERMICOMMON
