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

SUBROUTINE GENDIAG

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE DIAGARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: INFO
  INTEGER :: I, J, K, M, LWORK
  REAL(LATTEPREC), ALLOCATABLE :: GENWORK(:), STMP(:,:)
  IF (EXISTERROR) RETURN

  LWORK = 3*HDIM - 1
  ALLOCATE(GENWORK(LWORK), STMP(HDIM, HDIM))

  EVECS = H
  STMP = SMAT

  CALL DSYGV(1,'V','U', HDIM, EVECS, HDIM, STMP, HDIM, EVALS, GENWORK, &
       LWORK, INFO)

  IF (DEBUGON .EQ. 1) THEN

     DO I = 1, HDIM 
        PRINT*, I, EVALS(I)
     ENDDO

     DO I = 1, HDIM
        WRITE(6,'(100F12.6)') (EVECS(J,I), J = 1, HDIM)
     ENDDO

  ENDIF

  DEALLOCATE(GENWORK, STMP)

  RETURN

END SUBROUTINE GENDIAG
