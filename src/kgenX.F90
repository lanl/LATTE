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

SUBROUTINE KGENX

  USE CONSTANTS_MOD
  USE NONOARRAY
  USE KSPACEARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, INFO, II, LWORK,  LIWORK, LRWORK
  INTEGER, ALLOCATABLE :: IWORK(:)
  REAL(LATTEPREC) :: INVSQRT
  REAL(LATTEPREC), ALLOCATABLE :: EVAL(:), RWORK(:)
  COMPLEX(LATTEPREC), ALLOCATABLE :: UK(:,:), KTMPMAT(:,:), WORK(:)
  COMPLEX(LATTEPREC) :: ALPHA, BETA
  IF (EXISTERROR) RETURN


#ifdef XHEEV

  LWORK = 2*HDIM - 1
  LRWORK = 3*HDIM - 2

  ALLOCATE(WORK(LWORK), RWORK(3*HDIM - 2))

#elif defined(XHEEVD)

  LWORK = 2*HDIM + HDIM*HDIM
  LRWORK = 1 + 5*HDIM + 2*HDIM*HDIM
  LIWORK = 3 + 5*HDIM

  ALLOCATE(WORK(LWORK), RWORK(LRWORK), IWORK(LIWORK))

#else

  LWORK = 2*HDIM - 1
  LRWORK = 3*HDIM - 2

  ALLOCATE(WORK(LWORK), RWORK(3*HDIM - 2))
 
#endif  


  ALLOCATE(UK(HDIM, HDIM), EVAL(HDIM), KTMPMAT(HDIM, HDIM))
  !
  ! X = U s^-1/2 U^dag
  !

  ! Eigenvectors overwrite S (S = U)

  DO II = 1, NKTOT

     UK(:,:) = SK(:,:,II)


#ifdef XHEEV
!     PRINT*, "ZHEEV"

     CALL ZHEEV('V', 'U', HDIM, UK, HDIM, EVAL, WORK, LWORK, RWORK, INFO)

#elif defined(XHEEVD)

!     PRINT*, "ZHEEVD"
     CALL ZHEEVD('V', 'U', HDIM, UK, HDIM, EVAL, &
          WORK, LWORK, RWORK, LRWORK, IWORK, LIWORK, INFO)
     
     IF (INFO .NE. 0) CALL ZHEEV('V', 'U', HDIM, UK, HDIM, EVAL, WORK, LWORK, RWORK, INFO)

#else

     CALL ZHEEV('V', 'U', HDIM, UK, HDIM, EVAL, WORK, LWORK, RWORK, INFO)

#endif     

!     CALL ZHEEV('V', 'U', HDIM, UK, HDIM, EVAL, WORK, LWORK, RWORK, INFO)

     IF (EVAL(1) .LT. ZERO) THEN
        CALL ERRORS("kgenX","Eigenvalue of complex S matrix < 0: STOP!")
     ENDIF

     DO I = 1, HDIM

        INVSQRT = ONE/SQRT(EVAL(I))

        DO J = 1, HDIM
           KTMPMAT(J,I) = UK(J,I) * INVSQRT
        ENDDO

     ENDDO

     ALPHA = CMPLX(ONE, ZERO)
     BETA = CMPLX(ZERO, ZERO)

     CALL ZGEMM('N', 'C', HDIM, HDIM, HDIM, ALPHA, KTMPMAT, HDIM, UK, &
          HDIM, BETA, KXMAT(:,:,II), HDIM)

  ENDDO

  DEALLOCATE(UK, EVAL, KTMPMAT)
  IF (ALLOCATED(WORK)) DEALLOCATE(WORK)
  IF (ALLOCATED(IWORK)) DEALLOCATE(IWORK)
  IF (ALLOCATED(RWORK)) DEALLOCATE(RWORK)

  RETURN

END SUBROUTINE KGENX

