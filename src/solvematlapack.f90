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

SUBROUTINE SOLVEMATLAPACK

  USE CONSTANTS_MOD
  USE FERMICOMMON
  USE SETUPARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I
  INTEGER :: INFO
  REAL(LATTEPREC), ALLOCATABLE :: IPIV(:)
  IF (EXISTERROR) RETURN

  ALLOCATE(IPIV(HDIM))

  IF (SPINON .EQ. 0) THEN

#ifdef DOUBLEPREC
     CALL  DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          BO, HDIM, BO, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL  SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          BO, HDIM, BO, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - BO)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

     BO = X2

#ifdef DOUBLEPREC
     CALL DGESV(HDIM,  HDIM, A, HDIM, IPIV, BO, HDIM, INFO)
#elif defined(SINGLEPREC)
     CALL SGESV(HDIM,  HDIM, A, HDIM, IPIV, BO, HDIM, INFO)
#endif

  ELSE

#ifdef DOUBLEPREC
     CALL  DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHOUP, HDIM, RHOUP, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL  SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          RHOUP, HDIM, RHOUP, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - RHOUP)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

     RHOUP = X2

#ifdef DOUBLEPREC
     CALL DGESV(HDIM,  HDIM, A, HDIM, IPIV, RHOUP, HDIM, INFO)
#elif defined(SINGLEPREC)
     CALL SGESV(HDIM,  HDIM, A, HDIM, IPIV, RHOUP, HDIM, INFO)
#endif


#ifdef DOUBLEPREC
     CALL  DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHODOWN, HDIM, RHODOWN, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
     CALL  SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
          RHODOWN, HDIM, RHODOWN, HDIM, 0.0, X2, HDIM)
#endif

     A = TWO*(X2 - RHODOWN)

     DO I = 1, HDIM
        A(I,I) = A(I,I) + ONE
     ENDDO

     RHODOWN = X2

#ifdef DOUBLEPREC
     CALL DGESV(HDIM,  HDIM, A, HDIM, IPIV, RHODOWN, HDIM, INFO)
#elif defined(SINGLEPREC)
     CALL SGESV(HDIM,  HDIM, A, HDIM, IPIV, RHODOWN, HDIM, INFO)
#endif

  ENDIF

  DEALLOCATE(IPIV)

  RETURN

END SUBROUTINE SOLVEMATLAPACK





