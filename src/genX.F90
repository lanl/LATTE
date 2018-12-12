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

SUBROUTINE GENX

  USE CONSTANTS_MOD
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, INFO
  REAL(LATTEPREC) :: INVSQRT
  REAL(LATTEPREC), ALLOCATABLE :: IDENTITY(:,:)

  IF (EXISTERROR) RETURN

  !
  ! X = U s^-1/2 U^dag
  !

  ! Eigenvectors overwrite S (S = U)

  UMAT = SMAT

#ifdef XSYEV

#ifdef DOUBLEPREC
  CALL DSYEV("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
       NONO_LWORK,  INFO)
#elif defined(SINGLEPREC)
  CALL SSYEV("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
       NONO_LWORK, INFO)
#endif

#elif defined(XSYEVD)

#ifdef DOUBLEPREC
  CALL DSYEVD("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
       NONO_LWORK, NONO_IWORK, NONO_LIWORK, INFO)

  IF (INFO .NE. 0) THEN
       
     INFO = 0

     CALL DSYEV("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
          NONO_LWORK,  INFO)

  ENDIF


#elif defined(SINGLEPREC)
  CALL SSYEVD("V", "U", HDIM, UMAT, HDIM, NONO_EVALS, NONO_WORK, &
       NONO_LWORK, NONO_IWORK, NONO_LIWORK, INFO)
#endif

#endif


  IF ( NONO_EVALS(1) .LE. ZERO ) THEN
     CALL FITTINGOUTPUT(1)
     CALL PANIC
     CALL ERRORS("genX","Eigenvalues of overlap matrix <= 0")
     IF (EXISTERROR) RETURN 
  ENDIF

  !  PRINT*,  MINVAL(NONO_EVALS), MAXVAL(NONO_EVALS)

  DO I = 1, HDIM

     INVSQRT = ONE/SQRT(NONO_EVALS(I))

     DO J = 1, HDIM
        NONOTMP(J,I) = UMAT(J,I) * INVSQRT
     ENDDO

  ENDDO

#ifdef DOUBLEPREC
  CALL DGEMM('N', 'T', HDIM, HDIM, HDIM, ONE, &
       NONOTMP, HDIM, UMAT, HDIM, ZERO, XMAT, HDIM)
#elif defined(SINGLEPREC)
  CALL SGEMM('N', 'T', HDIM, HDIM, HDIM, ONE, &
       NONOTMP, HDIM, UMAT, HDIM, ZERO, XMAT, HDIM)
#endif

  IF (DEBUGON .EQ. 1) THEN

     ALLOCATE(IDENTITY(HDIM, HDIM))

     PRINT*, "Caution - you're writing to file the X matrix!"

     OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myX.dat")

     DO I = 1, HDIM
        WRITE(31,10) (XMAT(I,J), J = 1, HDIM)
     ENDDO

     CLOSE(31)

10   FORMAT(100G18.8)

     ! Let's also check the inverse

     !     CALL DGEMM( 'N', 'N', HDIM, HDIM, HDIM, ONE, &
     !             XMAT, HDIM, XMAT, HDIM, ZERO, XSQ, HDIM)
     !     CALL DGEMM( 'N', 'N', HDIM, HDIM, HDIM, ONE, &
     !             XSQ, HDIM, IDENTITY, HDIM, ZERO, NONOTMP, HDIM)

#ifdef DOUBLEPREC

     CALL DGEMM( 'T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, SMAT, HDIM, ZERO, NONOTMP, HDIM)
     CALL DGEMM( 'N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, IDENTITY, HDIM)

#elif defined(SINGLEPREC)

     CALL SGEMM( 'T', 'N', HDIM, HDIM, HDIM, ONE, &
          XMAT, HDIM, SMAT, HDIM, ZERO, NONOTMP, HDIM)
     CALL SGEMM( 'N', 'N', HDIM, HDIM, HDIM, ONE, &
          NONOTMP, HDIM, XMAT, HDIM, ZERO, IDENTITY, HDIM)

#endif

     OPEN(UNIT=31, STATUS="UNKNOWN", FILE="myXcheck.dat")

     DO I = 1, HDIM
        WRITE(31,10) (IDENTITY(I,J), J = 1, HDIM)
     ENDDO

     CLOSE(31)

     DEALLOCATE (IDENTITY)

  ENDIF


  RETURN

END SUBROUTINE GENX
