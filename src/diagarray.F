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

MODULE DIAGARRAY

  USE MYPRECISION 

  IMPLICIT NONE
  SAVE

  INTEGER :: DIAG_LWORK, DIAG_LIWORK, DIAG_LZWORK, DIAG_LRWORK
  INTEGER :: ZHEEVD_LWORK, ZHEEVD_LRWORK, ZHEEVD_LIWORK
  INTEGER, ALLOCATABLE :: DIAG_IWORK(:), IFAIL(:)
  INTEGER, ALLOCATABLE :: ZHEEVD_IWORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: DIAG_WORK(:), DIAG_RWORK(:)
  REAL(LATTEPREC), ALLOCATABLE ::  ZHEEVD_RWORK(:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: DIAG_ZWORK(:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: ZHEEVD_WORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: EVALS(:), EVECS(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: UPEVALS(:), UPEVECS(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: DOWNEVALS(:), DOWNEVECS(:,:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: ZBO(:,:), KHTMP(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: KEVALS(:,:), CPLIST(:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: KEVECS(:,:,:)
  REAL(LATTEPREC), PARAMETER :: EXPTOL = 30.0
  REAL(LATTEPREC) :: NUMLIMIT

  !  NUMLIMIT = EXP(-EXPTOL)

END MODULE DIAGARRAY
