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

SUBROUTINE NORMS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, LWORK, INFO
  REAL(LATTEPREC), ALLOCATABLE :: HP(:,:), X2(:,:), DBRHO(:,:), DBH(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: XTX(:,:), EVALS(:), WORK(:)
  REAL(LATTEPREC), ALLOCATABLE :: TRACEARR(:)
  REAL(LATTEPREC) :: FROBCOM, FROBIDEM, TRERR, TRTMP
  REAL(LATTEPREC) :: TWOCOM, TWOIDEM, TRRHOS


  LWORK = 3*HDIM - 1

  ALLOCATE(DBRHO(HDIM, HDIM), DBH(HDIM, HDIM))
  ALLOCATE(EVALS(HDIM), WORK(LWORK))

  IF (BASISTYPE .EQ. "ORTHO") THEN

     IF (LATTEPREC .EQ. 8) THEN

        DBRHO = BO/2.0D0
        DBH = H

     ELSEIF (LATTEPREC .EQ. 4) THEN

        DBRHO = DBLE(BO)/2.0D0
        DBH = DBLE(H)

     ENDIF


     ALLOCATE(HP(HDIM, HDIM))

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          DBH, HDIM, DBRHO, HDIM, 0.0D0, HP, HDIM)

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, -1.0D0, &
          DBRHO, HDIM, DBH, HDIM, 1.0D0, HP, HDIM)

     FROBCOM = 0.0D0

     DO I = 1, HDIM
        DO J = 1, HDIM

           FROBCOM = FROBCOM + HP(J,I)*HP(J,I)

        ENDDO
     ENDDO

     FROBCOM = SQRT(FROBCOM)

     ALLOCATE(XTX(HDIM, HDIM)) 

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          HP, HDIM, HP, HDIM, 0.0D0, XTX, HDIM)   

     CALL DSYEV('N', 'U', HDIM, XTX, HDIM, EVALS, WORK, LWORK, INFO)

     TWOCOM = SQRT(MAXVAL(EVALS))  

     DEALLOCATE(HP)

     ALLOCATE(X2(HDIM, HDIM))

     X2 = DBRHO

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          DBRHO, HDIM, DBRHO, HDIM, -1.0D0, X2, HDIM) 

     FROBIDEM = ZERO

     DO I = 1, HDIM
        DO J = 1, HDIM

           FROBIDEM = FROBIDEM + X2(J,I)*X2(J,I)

        ENDDO
     ENDDO

     FROBIDEM = SQRT(FROBIDEM)

     ! Here X2 = P^2 - P

     ! For the 2-norm we need the largest eigenvalue of X2^T X2

     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          X2, HDIM, X2, HDIM, 0.0D0, XTX, HDIM)   

     CALL DSYEV('N', 'U', HDIM, XTX, HDIM, EVALS, WORK, LWORK, INFO)

     TWOIDEM = SQRT(MAXVAL(EVALS))

     DEALLOCATE(X2, XTX)


     !  ALLOCATE(TRACEARR(HDIM))
     !  DO I = 1, HDIM
     !     TRACEARR(I) = DBRHO(I,I)
     !  ENDDO

     !  TRTMP = SUM(TRACEARR)

     TRTMP = 0.0D0

     DO I = 1, HDIM

        TRTMP = TRTMP + DBRHO(I,I)

     ENDDO

     TRERR = ABS(TRTMP - DBLE(BNDFIL)*DBLE(HDIM))




     !     DEALLOCATE(DBRHO, DBH, EVALS, WORK)
     WRITE(6,'("")')    
     WRITE(6,'("HDIM = ", I6)') HDIM
     WRITE(6,'("Occupation error  = ", G26.16)') TRERR
     WRITE(6,'("Frobenius norms:")')
     WRITE(6,'("Commutation error = ", G26.16)') FROBCOM
     WRITE(6,'("Idempotency error = ", G26.16)') FROBIDEM
     WRITE(6,'("Two norms:")')
     WRITE(6,'("Commutation error = ", G26.16)') TWOCOM
     WRITE(6,'("Idempotency error = ", G26.16)') TWOIDEM
     WRITE(6,10) HDIM, TRERR, FROBCOM, FROBIDEM, TWOCOM, TWOIDEM

10   FORMAT(I8, 5G26.16)

  ELSEIF (BASISTYPE .EQ. "NONORTHO") THEN

     ! Check rho S rho:


     IF (LATTEPREC .EQ. 8) THEN

        DBRHO = BO/2.0D0
        DBH = H

     ELSEIF (LATTEPREC .EQ. 4) THEN

        DBRHO = DBLE(BO)/2.0D0
        DBH = DBLE(H)

     ENDIF

     ALLOCATE(XTX(HDIM, HDIM))
     ALLOCATE(X2(HDIM, HDIM))

     !     X2 = DBRHO

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          DBRHO, HDIM, SMAT, HDIM, 0.0D0, XTX, HDIM)

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          XTX, HDIM,  DBRHO, HDIM, 0.0D0, X2, HDIM)

     X2 = X2 - DBRHO


     CALL DGEMM('T', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          X2, HDIM, X2, HDIM, 0.0D0, XTX, HDIM)   

     CALL DSYEV('N', 'U', HDIM, XTX, HDIM, EVALS, WORK, LWORK, INFO)

     TWOIDEM = SQRT(MAXVAL(EVALS))

     TRRHOS = 0.0D0

     DO I = 1, HDIM
        DO J = 1, HDIM
           TRRHOS = TRRHOS + DBRHO(J,I)*SMAT(I,J)
        ENDDO
     ENDDO

     PRINT*, "Idempotency error: || rSr - r ||_2 =", TWOIDEM
     PRINT*, "Tr(rho S) =", TRRHOS

     DEALLOCATE(X2, XTX)
  ENDIF

  DEALLOCATE(DBRHO, DBH, EVALS, WORK)
  RETURN

END SUBROUTINE NORMS
