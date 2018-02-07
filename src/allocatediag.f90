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

SUBROUTINE ALLOCATEDIAG

  USE CONSTANTS_MOD
  USE DIAGARRAY
  USE KSPACEARRAY

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  IF (KON .EQ. 0) THEN


#ifdef XSYEV

     DIAG_LWORK = 3*HDIM - 1
     ALLOCATE(DIAG_WORK(DIAG_LWORK))

#elif defined(XSYEVD)

     DIAG_LWORK = 1 + 6*HDIM + 2*HDIM*HDIM
     DIAG_LIWORK = 3 + 5*HDIM

     ALLOCATE(DIAG_WORK(DIAG_LWORK), DIAG_IWORK(DIAG_LIWORK))

#endif


     IF (SPINON .EQ. 0) THEN
        ALLOCATE (EVALS(HDIM), EVECS(HDIM, HDIM))
     ELSE
        ALLOCATE (UPEVALS(HDIM), UPEVECS(HDIM, HDIM))
        ALLOCATE (DOWNEVALS(HDIM), DOWNEVECS(HDIM, HDIM))
     ENDIF

  ELSE ! Allocate arrays if we're doing k-space

     ! Using ZHEEV..


     DIAG_LZWORK = 2*HDIM-1
     DIAG_LRWORK = 3*HDIM-2

     ALLOCATE(DIAG_RWORK(DIAG_LRWORK), DIAG_ZWORK(DIAG_LZWORK))

     ! Using ZHEEVD

     !     ZHEEVD_LWORK = 2*HDIM + HDIM*HDIM
     !     ZHEEVD_LRWORK = 1 + 5*HDIM + 2*HDIM*HDIM
     !     ZHEEVD_LIWORK = 3 + 5*HDIM

     !     ALLOCATE(ZHEEVD_WORK(ZHEEVD_LWORK), ZHEEVD_RWORK(ZHEEVD_LRWORK), &
     !              ZHEEVD_IWORK(ZHEEVD_LIWORK))


     !     ALLOCATE(KHTMP(HDIM, HDIM))

     IF (SPINON .EQ. 0) THEN
        ALLOCATE (KEVALS(HDIM, NKTOT), KEVECS(HDIM, HDIM, NKTOT), &
             CPLIST(HDIM*NKTOT))
     ENDIF

  ENDIF


  NUMLIMIT = EXP(-EXPTOL)

  RETURN

END SUBROUTINE ALLOCATEDIAG
