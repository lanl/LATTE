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

SUBROUTINE SP2FERMI

  !
  ! This is an implementation of Niklasson's SP2 algorithm for the
  ! Fermi operator (i.e., finite temperature SP2)
  !
  ! GERSHGORIN and SP2FERMIINIT must be run first to initialize
  !

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

#ifdef GPUON

  CALL sp2fermi_gpu(bndfil, hdim, spinon, bo, rhoup, rhodown, maxeval, &
       h, hup, hdown, maxminusmin, chempot, norecs, signlist, beta0, &
       breaktol, latteprec)

#elif defined(GPUOFF)


  INTEGER :: I, J, II
  INTEGER :: ITER, BREAKLOOP
  REAL(LATTEPREC) :: OCC, GERSHFACT
  REAL(LATTEPREC) :: TRX, TRX1, TRXOMX, TRX2
  REAL(LATTEPREC) :: LAMBDA
  REAL(LATTEPREC) :: MAXSHIFT = ONE
  REAL(LATTEPREC) :: OCCERROR
  REAL(LATTEPREC) :: PREVERROR, PREVERROR2, PREVERROR3
  REAL(LATTEPREC) :: GEMM_ALPHA, GEMM_BETA
  IF (EXISTERROR) RETURN

  ITER = 0

  BREAKLOOP = 0

  PREVERROR = ZERO
  PREVERROR2 = ZERO
  PREVERROR3 = ZERO

  IF (SPINON .EQ. 0) THEN

     OCC = BNDFIL*FLOAT(HDIM)

     DO WHILE (BREAKLOOP .EQ. 0 .AND. ITER .LT. 50)

        ITER = ITER + 1

        BO = -H/MAXMINUSMIN

        GERSHFACT = (MAXEVAL - CHEMPOT)/MAXMINUSMIN

        DO I = 1, HDIM

           BO(I,I) = GERSHFACT + BO(I,I)

        ENDDO


        DO II = 1, NORECS

           X2 = BO

           IF (SIGNLIST(II) .EQ. 1) THEN

#ifdef DOUBLEPREC
              CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, -1.0D0, &
                   X2, HDIM, X2, HDIM, 2.0D0, BO, HDIM)
#elif defined(SINGLEPREC)
              CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, -1.0, &
                   X2, HDIM, X2, HDIM, 2.0, BO, HDIM)
#endif

           ELSE

#ifdef DOUBLEPREC
              CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
                   X2, HDIM, X2, HDIM, 0.0D0, BO, HDIM)
#elif defined(SINGLEPREC)
              CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
                   X2, HDIM, X2, HDIM, 0.0, BO, HDIM)
#endif

           ENDIF

        ENDDO

        TRX = ZERO
        TRX2 = ZERO

        !
        ! Now we're getting TrX1 where X1 = beta0*(X-X^2)
        !

        DO I = 1, HDIM

           TRX = TRX + BO(I,I)

           DO J = 1, HDIM

              TRX2 = TRX2 + BO(J,I)*BO(J,I)

           ENDDO

        ENDDO



        TRXOMX = TRX - TRX2

	!PRINT*, TRX, TRX2

        LAMBDA = (OCC - TRX)/(BETA0*TRXOMX)

        !
        ! New CHEMPOT
        !

        IF (ABS(LAMBDA) .GT. MAXSHIFT) LAMBDA = SIGN(MAXSHIFT, LAMBDA)

        CHEMPOT = CHEMPOT + LAMBDA

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR
        OCCERROR = ABS(OCC - TRX)

#ifdef DOUBLEPREC

        IF (OCCERROR .LT. BREAKTOL) THEN
           BREAKLOOP = 1
        ENDIF

#elif defined(SINGLEPREC)

        IF (OCCERROR .EQ. PREVERROR .OR. &
             OCCERROR .EQ. PREVERROR2 .OR. &
             OCCERROR .EQ. PREVERROR3 .OR. ITER .EQ. 10 ) THEN

           BREAKLOOP = 1

        ENDIF

#endif

     ENDDO

     IF (ITER .EQ. 100) THEN
        CALL ERRORS("sp2fermi","SP2FERMI is not converging: STOP!")
     ENDIF

     !
     ! If you forget the following you'll spend about a day
     ! trying to find the bug in every other subroutine...
     !

     !     BO = TWO*BO

     ! Update occupancy

     ! bo = 2 x ( bo + lambda*(beta0*(bo * (I - bo)))) <- we put this into
     ! GEMM-form

     GEMM_ALPHA = -TWO*LAMBDA*BETA0
     GEMM_BETA = TWO*(ONE + LAMBDA*BETA0)

     X2 = BO

#ifdef DOUBLEPREC
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, GEMM_ALPHA, &
          X2, HDIM, X2, HDIM, GEMM_BETA, BO, HDIM)
#elif defined(SINGLEPREC)
     CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, GEMM_ALPHA, &
          X2, HDIM, X2, HDIM, GEMM_BETA, BO, HDIM)
#endif

  ELSEIF (SPINON .EQ. 1) THEN

     DO WHILE ( BREAKLOOP .EQ. 0 )

        ITER = ITER  + 1

        IF (ITER .EQ. 50) THEN
           CALL ERRORS("sp2fermi","SP2FERMI is not converging: STOP")
        ENDIF

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 RHOUP(I,I) = (MAXEVAL - HUP(I,I) - CHEMPOT)/MAXMINUSMIN

                 RHODOWN(I,I) = (MAXEVAL - HDOWN(I,I) - CHEMPOT) / &
                      MAXMINUSMIN

              ELSE

                 RHOUP(J,I) = (ZERO - HUP(J,I))/MAXMINUSMIN
                 RHOUP(I,J) = RHOUP(J,I)

                 RHODOWN(J,I) = (ZERO - HDOWN(J,I))/MAXMINUSMIN
                 RHODOWN(I,J) = RHODOWN(J,I)

              ENDIF

           ENDDO
        ENDDO


        DO II = 1, NORECS

           !
           ! X*X
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

           RHOUP = RHOUP + SIGNLIST(II)*(RHOUP - X2UP)
           RHODOWN = RHODOWN + SIGNLIST(II)*(RHODOWN - X2DOWN)

        ENDDO

        TRX = ZERO
        TRX1 = ZERO

        !
        ! Now we're getting TrX1 where X1 = beta0*(X-X^2)
        !

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 TRX1 = TRX1 + RHOUP(I,I)*(ONE - RHOUP(I,I)) + &
                      RHODOWN(I,I)*(ONE - RHODOWN(I,I))

              ELSE

                 TRX1 = TRX1 - TWO*(RHOUP(J,I)*RHOUP(J,I) + &
                      RHODOWN(J,I)*RHODOWN(J,I))

              ENDIF

           ENDDO

           TRX = TRX + RHOUP(I,I) + RHODOWN(I,I)

        ENDDO

        TRX1 = BETA0*TRX1

        LAMBDA = (TOTNE - TRX)/TRX1

        !
        ! New CHEMPOT
        !

        IF (ABS(LAMBDA) .GT. MAXSHIFT) THEN
           LAMBDA = SIGN(MAXSHIFT, LAMBDA)
        ENDIF

        CHEMPOT = CHEMPOT + LAMBDA

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR
        OCCERROR = ABS(TOTNE - TRX)

        !
        ! Convergence tests for double and single precision runs
        !

        !        PRINT*, ITER, OCCERROR, PREVERROR, PREVERROR2, PREVERROR3

#ifdef DOUBLEPREC

        IF (OCCERROR .LT. BREAKTOL) THEN
           BREAKLOOP = 1
        ENDIF

#elif defined(SINGLEPREC)

        IF ( (ITER .GE. 2) .AND. (OCCERROR .EQ. PREVERROR .OR. &
             OCCERROR .EQ. PREVERROR2 .OR. &
             OCCERROR .EQ. PREVERROR3 .OR. &
             ITER .EQ. 10) ) THEN

           BREAKLOOP = 1

        ENDIF

#endif

     ENDDO

  ENDIF

#endif

  RETURN

END SUBROUTINE SP2FERMI


!  NEWBETA = (TRXPLUS - TRXMINUS)/(TWO*DELTA*TRXOMX)

!  PRINT*, "SP2FERMI: BETA = ", NEWBETA
!  PRINT*, "SP2FERMI: KBT = ", ONE/NEWBETA
