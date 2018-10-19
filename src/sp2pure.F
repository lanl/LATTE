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

SUBROUTINE SP2PURE

  !
  ! This subroutine implements Niklasson's SP2 density matrix purification
  ! algorithm.
  !


  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPINARRAY
  USE NONOARRAY
  USE MYPRECISION

  IMPLICIT NONE

#ifdef GPUON

  INTEGER :: sp2convint

  IF (SP2CONV .EQ. "REL") THEN
     sp2convint = 0
  ELSEIF (SP2CONV .EQ. "ABS") THEN
     sp2convint = 1
  ELSE
     CALL ERRORS("sp2pure",'("Convergence criterion for SP2 not defined")')
  ENDIF

  IF (BASISTYPE .EQ. "ORTHO") THEN

     CALL sp2purify(bndfil, hdim, spinon, bo, rhoup, rhodown, maxeval, &
          h, hup, hdown, maxminusmin, minsp2iter, sp2convint, latteprec)

  ELSE

     CALL sp2purify(bndfil, hdim, spinon, bo, rhoup, rhodown, maxeval, &
          orthoh, orthohup, orthohdown, maxminusmin, minsp2iter, &
          sp2convint, latteprec)

  ENDIF

  ! call sp2evolution(bndfil, hdim, spinon, bo, rhoup, rhodown, maxeval, &
  ! h, hup, hdown, maxminusmin, minsp2iter, sp2convint, latteprec)

#elif defined(GPUOFF)

  INTEGER :: I, J, ITER
  INTEGER :: BREAKLOOP
  INTEGER :: PUR
  REAL(LATTEPREC) :: TRX, OCC, TRX2, LIMDIFF
  REAL(LATTEPREC) :: GERSHFACT
  REAL(LATTEPREC) :: IDEMPERR, TRXOLD
  REAL(LATTEPREC) :: IDEMPERR1, IDEMPERR2
#ifdef DOUBLEPREC
  REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0D-14
#elif defined(SINGLEPREC)
  REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0E-14
  IF (EXISTERROR) RETURN
#endif

  IDEMPERR = ZERO
  IDEMPERR1 = ZERO
  IDEMPERR2 = ZERO

  !
  ! We're also using Niklasson's scheme to determine convergence
  !

  OCC = BNDFIL*FLOAT(HDIM)


  IF (SPINON .EQ. 0) THEN

     ! Using intrinsics is probably better than coding this ourselves

     ! Build the starting guess

     IF (BASISTYPE .EQ. "ORTHO") THEN
        BO = -H/MAXMINUSMIN
     ELSE
        BO = -ORTHOH/MAXMINUSMIN
     ENDIF

     GERSHFACT =  MAXEVAL/MAXMINUSMIN

     TRX = ZERO

     DO I = 1, HDIM

        BO(I,I) = GERSHFACT + BO(I,I)
	TRX = TRX + BO(I,I)

     ENDDO



     ITER = 0

     BREAKLOOP = 0

     DO WHILE ( BREAKLOOP .EQ. 0 .AND. ITER .LT. 100 )

        ITER = ITER + 1

        X2 = BO

#ifdef DOUBLEPREC

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             BO, HDIM, BO, HDIM, ONE, X2, HDIM)

#elif defined(SINGLEPREC)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             BO, HDIM, BO, HDIM, ONE, X2, HDIM)

#endif

        TRX2 = ZERO

        DO I = 1, HDIM
           TRX2 = TRX2 + X2(I,I)
        ENDDO


	LIMDIFF = ABS(TRX - TRX2 - OCC) - ABS(TRX + TRX2 - OCC)

        IF ( LIMDIFF .GE. IDEMTOL ) THEN

	   BO = BO + X2

           TRX = TRX + TRX2

        ELSEIF ( LIMDIFF .LT. -IDEMTOL ) THEN

	   BO = BO - X2

           TRX = TRX - TRX2

        ENDIF

	IDEMPERR2 = IDEMPERR1
	IDEMPERR1 = IDEMPERR
        IDEMPERR = ABS(TRX2)

        !        PRINT*, ITER, IDEMPERR, ABS(TRX - OCC)
        !	WRITE(*,10) ITER, IDEMPERR, IDEMPERR2 - IDEMPERR
10      FORMAT(I4, 2G30.18)
 	IF (SP2CONV .EQ. "REL" .AND. ITER .GE. MINSP2ITER &
             .AND. (IDEMPERR2 .LE. IDEMPERR .OR. &
             IDEMPERR .LT. IDEMTOL)) BREAKLOOP = 1

  !	IF (ITER .EQ. 30) BREAKLOOP=1
        IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LT. IDEMTOL) BREAKLOOP = 1

     ENDDO

     !     PRINT*, "TRX = ", TRX

     IF (ITER .EQ. 100) THEN
        CALL PANIC
        CALL ERRORS("sp2pure","SP2 purification is not converging: STOP!")
     ENDIF

     BO = TWO*BO

  ELSE

     !
     ! Start by remapping the two H matrices such that
     ! all eigenvalues are [0:1]. We have the Gersgorin bounds
     ! for both Hup and Hdown
     !

     IF (BASISTYPE .EQ. "ORTHO") THEN
        RHOUP = -HUP/MAXMINUSMIN
        RHODOWN = -HDOWN/MAXMINUSMIN
     ELSE
        RHOUP = -ORTHOHUP/MAXMINUSMIN
        RHODOWN = -ORTHOHDOWN/MAXMINUSMIN
     ENDIF

     GERSHFACT = MAXEVAL/MAXMINUSMIN

     TRX = ZERO
     DO I = 1, HDIM
        RHOUP(I,I) = GERSHFACT + RHOUP(I,I)
        RHODOWN(I,I) = GERSHFACT + RHODOWN(I,I)
        TRX = TRX + RHOUP(I,I) + RHODOWN(I,I)
     ENDDO

     !     PRINT*, "TRX = ", TRX, "TOTNE = ", TOTNE

     ITER = 0

     BREAKLOOP = 0

     DO WHILE ( BREAKLOOP .EQ. 0 .AND. ITER .LT. 100 )

        ITER = ITER + 1

        X2UP = RHOUP
        X2DOWN = RHODOWN

#ifdef DOUBLEPREC

        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             RHOUP, HDIM, RHOUP, HDIM, ONE, X2UP, HDIM)
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             RHODOWN, HDIM, RHODOWN, HDIM, ONE, X2DOWN, HDIM)

#elif defined(SINGLEPREC)

        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             RHOUP, HDIM, RHOUP, HDIM, ONE, X2UP, HDIM)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
             RHODOWN, HDIM, RHODOWN, HDIM, ONE, X2DOWN, HDIM)

#endif

        TRX2 = ZERO

        DO I = 1, HDIM
           TRX2 = TRX2 + X2UP(I,I) + X2DOWN(I,I)
        ENDDO

	TRXOLD = TRX

	LIMDIFF = ABS(TRX - TRX2 - TOTNE) - ABS(TRX + TRX2 - TOTNE)

        IF (LIMDIFF .GE. IDEMTOL) THEN

           !$OMP PARALLEL DO DEFAULT(NONE) &
           !$OMP SHARED(RHOUP, RHODOWN, X2UP, X2DOWN, HDIM) &
           !$OMP PRIVATE(I,J)

           DO I = 1, HDIM
              DO J = 1, HDIM

                 RHOUP(J,I) = RHOUP(J,I) + X2UP(J,I)
                 RHODOWN(J,I) = RHODOWN(J,I) + X2DOWN(J,I)

              ENDDO
           ENDDO

           !$OMP END PARALLEL DO
           !$OMP BARRIER

           TRX = TRX + TRX2

        ELSEIF ( LIMDIFF .LT. -IDEMTOL ) THEN


           !$OMP PARALLEL DO DEFAULT(NONE) &
           !$OMP SHARED(RHOUP, RHODOWN, X2UP, X2DOWN, HDIM) &
           !$OMP PRIVATE(I,J)

           DO I = 1, HDIM
              DO J = 1, HDIM

                 RHOUP(J,I) = RHOUP(J,I) - X2UP(J,I)
                 RHODOWN(J,I) = RHODOWN(J,I) - X2DOWN(J,I)

              ENDDO
           ENDDO

           !$OMP END PARALLEL DO
           !$OMP BARRIER

           TRX = TRX - TRX2

        ENDIF

	IDEMPERR2 = IDEMPERR1
	IDEMPERR1 = IDEMPERR
        IDEMPERR = ABS(TRX2)

        !        PRINT*, ITER, IDEMPERR, LIMDIFF, TRX - TOTNE

 	IF (SP2CONV .EQ. "REL" .AND. ITER .GE. MINSP2ITER &
             .AND. (IDEMPERR2 .LE. IDEMPERR .OR. &
	     IDEMPERR .LT. IDEMTOL) ) BREAKLOOP = 1

        IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LE. IDEMTOL) BREAKLOOP = 1

     ENDDO

     IF (ITER .EQ. 100) THEN
        CALL PANIC
        CALL ERRORS("sp2pure","SP2 purification not converging: STOP!")
     ENDIF

  ENDIF

#endif

  RETURN

END SUBROUTINE SP2PURE
