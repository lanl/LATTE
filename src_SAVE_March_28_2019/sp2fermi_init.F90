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

SUBROUTINE SP2FERMIINIT

  !
  ! This is an implementation of Niklasson's SP2 algorithm for the
  ! Fermi operator (i.e., finite temperature SP2)
  !

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, II,  ITER
  INTEGER :: BREAKLOOP, MYSIGN
  REAL(LATTEPREC) :: TRX, OCC
  REAL(LATTEPREC) :: TRX2, TRXOMX, TRX1
  REAL(LATTEPREC) :: LAMBDA
  REAL(LATTEPREC) :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: OCCERROR = ZERO
  REAL(LATTEPREC), ALLOCATABLE :: X0X1(:,:), X1X0(:,:), X1(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: X0X1UP(:,:), X0X1DOWN(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: X1X0UP(:,:), X1X0DOWN(:,:)
  REAL(LATTEPREC), ALLOCATABLE :: X1UP(:,:), X1DOWN(:,:)
  REAL(LATTEPREC) :: PREVERROR, PREVERROR2, PREVERROR3
  IF (EXISTERROR) RETURN
  !  REAL(LATTEPREC) :: DELTA = 0.01D0, TRXPLUS, TRXMINUS, NEWBETA, NEWCHEMPOT

  !
  ! We'll have spin and non-spin dependent versions separate
  !

  ITER = 0

  BREAKLOOP = 0

  PREVERROR = ZERO
  PREVERROR2 = ZERO
  PREVERROR3 = ZERO

  IF (SPINON .EQ. 0) THEN

     ALLOCATE( X0X1(HDIM, HDIM), X1X0(HDIM, HDIM), X1(HDIM, HDIM) )

     OCC = BNDFIL*FLOAT(HDIM)

     DO WHILE ( BREAKLOOP .EQ. 0 )

        ITER = ITER + 1

        IF (ITER .EQ. 100) THEN
           CALL ERRORS("sp2fermi_init","SP2FERMIINIT is not converging: STOP!")
        ENDIF

        X1 = ZERO

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 BO(I,I) = (MAXEVAL - H(I,I) - CHEMPOT)/MAXMINUSMIN

              ELSE

                 BO(J,I) = (ZERO - H(J,I))/MAXMINUSMIN
                 BO(I,J) = BO(J,I)

              ENDIF

           ENDDO

           X1(I,I) = MINUSONE/MAXMINUSMIN

        ENDDO

        DO II = 1, NORECS

           !
           ! Density matrix squared
           !

#ifdef DOUBLEPREC
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                BO, HDIM, BO, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                BO, HDIM, BO, HDIM, 0.0, X2, HDIM)
#endif

           TRX = ZERO
           TRX2 = ZERO
           DO I = 1, HDIM
              TRX = TRX + BO(I,I)
              TRX2 = TRX2 + X2(I,I)
           ENDDO

           IF ( ABS(TRX2 - OCC) .LT. ABS(TWO*TRX - TRX2 - OCC) ) THEN

              MYSIGN = -1

           ELSE

              MYSIGN = 1

           ENDIF

           SIGNLIST(II) = MYSIGN

           !
           ! Density matrix X reponse
           !

#ifdef DOUBLEPREC
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                BO, HDIM, X1, HDIM, 0.0D0, X0X1, HDIM)
#elif defined(SINGLEPREC)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                BO, HDIM, X1, HDIM, 0.0, X0X1, HDIM)
#endif

#ifdef DOUBLEPREC
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                X1, HDIM, BO, HDIM, 0.0D0, X1X0, HDIM)
#elif defined(SINGLEPREC)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                X1, HDIM, BO, HDIM, 0.0, X1X0, HDIM)
#endif

           !
           ! Update response
           !


           X1 = X1 + FLOAT(MYSIGN)*(X1 - X0X1 - X1X0)

           !
           ! Update density matrix
           !

           BO = BO + FLOAT(MYSIGN)*(BO - X2)

        ENDDO

        TRX = ZERO
        TRXOMX = ZERO
        TRX1 = ZERO

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 TRXOMX = TRXOMX + BO(I,I)*(ONE - BO(I,I))

              ELSE

                 TRXOMX = TRXOMX - TWO*BO(J,I)*BO(J,I)

              ENDIF

           ENDDO

           TRX1 = TRX1 + X1(I,I)
           TRX = TRX + BO(I,I)

        ENDDO

        BETA0 = TRX1/TRXOMX

        KBT = ABS(ONE/BETA0)

        LAMBDA = (OCC - TRX)/TRX1

        IF (ABS(LAMBDA) .GT. MAXSHIFT) THEN
           LAMBDA = SIGN(MAXSHIFT, LAMBDA)
        ENDIF

        CHEMPOT = CHEMPOT + LAMBDA

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR

        OCCERROR = ABS(OCC - TRX)

        ! How we figure if we've reached convergence. An absolute
        ! tolerance works well in double precision, but in single
        ! precision we need to look for noise when we're near
        ! self-consistency

#ifdef DOUBLEPREC

        IF (OCCERROR .LT. BREAKTOL) THEN

           BREAKLOOP = 1

        ENDIF

#elif defined(SINGLEPREC)

        IF (OCCERROR .EQ. PREVERROR .OR. &
             OCCERROR .EQ. PREVERROR2 .OR. &
             OCCERROR .EQ. PREVERROR3 .OR. ITER .EQ. 25 ) THEN

           BREAKLOOP = 1

        ENDIF

#endif

     ENDDO

     !
     ! Little bit of fine tuning...
     !

     BO = TWO*(BO + LAMBDA*X1)

     DEALLOCATE( X0X1, X1X0 , X1 )

  ELSE

     ALLOCATE(X0X1UP(HDIM, HDIM), X0X1DOWN(HDIM, HDIM))
     ALLOCATE(X1X0UP(HDIM, HDIM), X1X0DOWN(HDIM, HDIM))
     ALLOCATE(X1UP(HDIM, HDIM), X1DOWN(HDIM, HDIM))

     DO WHILE (BREAKLOOP .EQ. 0)

        ITER = ITER + 1

        IF (ITER .EQ. 100) THEN
           CALL ERRORS("sp2fermi_init","SP2FERMIINIT is not converging: STOP!")
        ENDIF

        X1UP = ZERO

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 RHOUP(I,I) = (MAXEVAL - HUP(I,I) - CHEMPOT)/MAXMINUSMIN
                 RHODOWN(I,I) = (MAXEVAL - HDOWN(I,I) - CHEMPOT)/MAXMINUSMIN

              ELSE

                 RHOUP(J,I) = (ZERO - HUP(J,I))/MAXMINUSMIN
                 RHOUP(I,J) = RHOUP(J,I)

                 RHODOWN(J,I) = (ZERO - HDOWN(J,I))/MAXMINUSMIN
                 RHODOWN(I,J) = RHODOWN(J,I)

              ENDIF

           ENDDO

           X1UP(I,I) = MINUSONE/MAXMINUSMIN

        ENDDO

        X1DOWN = X1UP

        DO II = 1, NORECS

           !
           ! Density matrix squared
           !

#ifdef DOUBLEPREC
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                RHOUP, HDIM, RHOUP, HDIM, 0.0D0, X2UP, HDIM)
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                RHODOWN, HDIM, RHODOWN, HDIM, 0.0D0, X2DOWN, HDIM)
#elif defined(SINGLEPREC)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                RHOUP, HDIM, RHOUP, HDIM, 0.0, X2UP, HDIM)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                RHODOWN, HDIM, RHODOWN, HDIM, 0.0, X2DOWN, HDIM)
#endif

           TRX = ZERO
           TRX2 = ZERO
           DO I = 1, HDIM
              TRX = TRX + RHOUP(I,I) + RHODOWN(I,I)
              TRX2 = TRX2 + X2UP(I,I) + X2DOWN(I,I)
           ENDDO

           IF (ABS(TRX2 - TOTNE) .LT. ABS(TWO*TRX - TRX2 - TOTNE)) THEN

              MYSIGN = -1

           ELSE

              MYSIGN = 1

           ENDIF

           SIGNLIST(II) = MYSIGN

#ifdef DOUBLEPREC

           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                RHOUP, HDIM, X1UP, HDIM, 0.0D0, X0X1UP, HDIM)
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                RHODOWN, HDIM, X1DOWN, HDIM, 0.0D0, X0X1DOWN, HDIM)

           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                X1UP, HDIM, RHOUP, HDIM, 0.0D0, X1X0UP, HDIM)
           CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                X1DOWN, HDIM, RHODOWN, HDIM, 0.0D0, X1X0DOWN, HDIM)

#elif defined(SINGLEPREC)

           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                RHOUP, HDIM, X1UP, HDIM, 0.0, X0X1UP, HDIM)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                RHODOWN, HDIM, X1DOWN, HDIM, 0.0, X0X1DOWN, HDIM)

           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                X1UP, HDIM, RHOUP, HDIM, 0.0, X1X0UP, HDIM)
           CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                X1DOWN, HDIM, RHODOWN, HDIM, 0.0, X1X0DOWN, HDIM)

#endif

           !
           ! Update response
           !


           X1UP = X1UP + MYSIGN*(X1UP - X0X1UP - X1X0UP)
           X1DOWN = X1DOWN + MYSIGN*(X1DOWN - X0X1DOWN - X1X0DOWN)

           !
           ! Update density matrix
           !

           RHOUP = RHOUP + MYSIGN*(RHOUP - X2UP)
           RHODOWN = RHODOWN + MYSIGN*(RHODOWN - X2DOWN)

        ENDDO

        TRX = ZERO
        TRXOMX = ZERO
        TRX1 = ZERO

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 TRXOMX = TRXOMX + RHOUP(I,I)*(ONE - RHOUP(I,I)) + &
                      RHODOWN(I,I)*(ONE - RHODOWN(I,I))

              ELSE

                 TRXOMX = TRXOMX - TWO*(RHOUP(J,I)*RHOUP(J,I) + &
                      RHODOWN(J,I)*RHODOWN(J,I))

              ENDIF

           ENDDO

           TRX1 = TRX1 + X1UP(I,I) + X1DOWN(I,I)
           TRX = TRX + RHOUP(I,I) + RHODOWN(I,I)

        ENDDO

        BETA0 = TRX1/TRXOMX

        KBT = ABS(ONE/BETA0)

        LAMBDA = (TOTNE - TRX)/TRX1

        IF (ABS(LAMBDA) .GT. MAXSHIFT) THEN
           LAMBDA = SIGN(MAXSHIFT, LAMBDA)
        ENDIF

        CHEMPOT = CHEMPOT + LAMBDA

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR

        OCCERROR = ABS(TOTNE - TRX)

        ! How we figure if we've reached convergence. An absolute
        ! tolerance works well in double precision, but in single
        ! precision we need to look for noise when we're near
        ! self-consistency

#ifdef DOUBLEPREC

        IF (OCCERROR .LT. BREAKTOL) THEN

           BREAKLOOP = 1

        ENDIF

#elif defined(SINGLEPREC)

        IF (OCCERROR .EQ. PREVERROR .OR. &
             OCCERROR .EQ. PREVERROR2 .OR. &
             OCCERROR .EQ. PREVERROR3 .OR. ITER .EQ. 25 ) THEN

           BREAKLOOP = 1

        ENDIF

#endif

     ENDDO

     !
     ! Little bit of fine tuning...
     !

     RHOUP = RHOUP + LAMBDA*X1UP
     RHODOWN = RHODOWN + LAMBDA*X1DOWN

     DEALLOCATE( X0X1UP, X0X1DOWN, X1X0UP, X1X0DOWN, X1UP, X1DOWN )

  ENDIF

  WRITE(6,'("# SP2FERMI: kbT in eV = ", F14.8)') KBT

  RETURN

END SUBROUTINE SP2FERMIINIT


!
! Now find the correct beta for calculating the corresponding entropy
!

!  ALLOCATE(X(HDIM, HDIM))

!  X = HALF*BO

!  TRXOMX = ZERO

!    DO I = 1, HDIM
!     DO J = I, HDIM

!        IF (I .EQ. J) THEN

!           TRXOMX = TRXOMX + X(I,I)*(ONE - X(I,I))

!        ELSE

!           TRXOMX = TRXOMX - TWO*X(J,I)*X(J,I)

!        ENDIF

!     ENDDO
!  ENDDO

!  NEWCHEMPOT = CHEMPOT + DELTA

!  DO I = 1, HDIM
!     DO J = I, HDIM

!        IF (I .EQ. J) THEN

!           X(I,I) = (MAXEVAL - H(I,I) - NEWCHEMPOT)/MAXMINUSMIN

!        ELSE

!           X(J,I) = (ZERO - H(J,I))/MAXMINUSMIN
!           X(I,J) = X(J,I)

!        ENDIF

!     ENDDO
!  ENDDO

!  DO II = 1, NORECS

!
! BO^2
!

!     IF (LATTEPREC .EQ. KIND(0.0D0)) THEN
!        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
!             X, HDIM, X, HDIM, 0.0D0, X2, HDIM)
!     ELSE
!        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
!             X, HDIM, X, HDIM, 0.0, X2, HDIM)
!     ENDIF

!
! We purify using the sequence of operations defined
! in SP2FERMIINIT
!

!     X = X + SIGNLIST(II)*(X - X2)

!  ENDDO

!  TRXPLUS = ZERO
!  DO I = 1, HDIM
!     TRXPLUS = TRXPLUS + X(I,I)
!  ENDDO

!  NEWCHEMPOT = CHEMPOT - DELTA

!  DO I = 1, HDIM
!     DO J = I, HDIM
!
!        IF (I .EQ. J) THEN

!           X(I,I) = (MAXEVAL - H(I,I) - NEWCHEMPOT)/MAXMINUSMIN

!        ELSE

!           X(J,I) = (ZERO - H(J,I))/MAXMINUSMIN
!           X(I,J) = X(J,I)

!        ENDIF

!     ENDDO
!  ENDDO

!  DO II = 1, NORECS

!
! BO^2
!

!     IF (LATTEPREC .EQ. KIND(0.0D0)) THEN
!        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
!             X, HDIM, X, HDIM, 0.0D0, X2, HDIM)
!     ELSE
!        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
!             X, HDIM, X, HDIM, 0.0, X2, HDIM)
!     ENDIF

!
! We purify using the sequence of operations defined
! in SP2FERMIINIT
!

!     X = X + SIGNLIST(II)*(X - X2)

!  ENDDO

!  TRXMINUS = ZERO
!  DO I = 1, HDIM
!     TRXMINUS = TRXMINUS + X(I,I)
!  ENDDO

!  DEALLOCATE(X)

!  NEWBETA = (TRXPLUS - TRXMINUS)/(TWO*DELTA*TRXOMX)

!  PRINT*, "SP2FERMI: BETA = ", NEWBETA
!  PRINT*, "SP2FERMI: KBT = ", ONE/NEWBETA
