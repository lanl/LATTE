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

SUBROUTINE FERMIEXPANS

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE FERMICOMMON
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER
  INTEGER :: BREAKLOOP
  REAL(LATTEPREC) :: OCC, OCCERROR
  REAL(LATTEPREC) :: PREVERROR, PREVERROR2, PREVERROR3
  REAL(LATTEPREC) :: BOVER2M, TRX, TRXOMX
  REAL(LATTEPREC) :: SHIFTCP, HFACT
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = ONE
  IF (EXISTERROR) RETURN

  OCC = BNDFIL*FLOAT(HDIM)

  ITER = 0

  BREAKLOOP = 0

  BOVER2M = ONE/(KBT*(TWO**(2 + FERMIM)))

  SCFS = 0

  PREVERROR = ZERO
  PREVERROR2 = ZERO
  PREVERROR3 = ZERO

  DO WHILE (BREAKLOOP .EQ. 0)

     SCFS = SCFS + 1

     ITER = ITER + 1

     IF (ITER .EQ. 100) THEN
        CALL PANIC
        CALL ERRORS("fermiexpans","Fermi expansion is not converging: STOP!")
     ENDIF

     IF (SPINON .EQ. 0) THEN

        BO = -BOVER2M*H

	HFACT = HALF + BOVER2M*CHEMPOT

	DO I = 1, HDIM
           BO(I,I) = HFACT + BO(I,I)
	ENDDO

     ELSE

        RHOUP = -BOVER2M*HUP
        RHODOWN = -BOVER2M*HDOWN

        HFACT = HALF + BOVER2M*CHEMPOT

        DO I = 1, HDIM
           RHOUP(I,I) = HFACT + RHOUP(I,I)
           RHODOWN(I,I) = HFACT + RHODOWN(I,I)
        ENDDO

     ENDIF
#ifdef GPUOFF

     DO I = 1, FERMIM

        IF (CGORLIB .EQ. 0) THEN

           ! Call GESV-based solver

           CALL SOLVEMATLAPACK

        ELSE

           ! Call the conjugate gradient solver on the CPU

           CALL SOLVEMATCG

        ENDIF

     ENDDO

#elif defined(GPUON)

     !
     ! This calls Sanville's CUDA routines on the GPU
     ! Now modified by MJC to send down FERMIM too
     !

     CALL SOLVE_MATRIX_CG(BO, RHOUP, RHODOWN, HDIM, CGTOL2, &
          SPINON, LATTEPREC, FERMIM)

#endif

     ! Modifying chemical potential

     TRX = ZERO
     TRXOMX = ZERO

     IF (SPINON .EQ. 0) THEN

        DO I = 1, HDIM
           DO J = I, HDIM

              IF (I .EQ. J) THEN

                 TRXOMX = TRXOMX + BO(I,I)*(ONE - BO(I,I))

              ELSE

                 TRXOMX = TRXOMX - TWO*BO(J,I)*BO(J,I)

              ENDIF

           ENDDO

           TRX = TRX + BO(I,I)

        ENDDO

        SHIFTCP = KBT*(OCC - TRX)/TRXOMX

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR
        OCCERROR = ABS(OCC - TRX)
        !        PRINT*, ITER, OCCERROR

     ELSE

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

           TRX = TRX + RHOUP(I,I) + RHODOWN(I,I)

        ENDDO

        SHIFTCP = KBT*(TOTNE - TRX)/TRXOMX

        PREVERROR3 = PREVERROR2
        PREVERROR2 = PREVERROR
        PREVERROR = OCCERROR
        OCCERROR = ABS(TOTNE - TRX)

     ENDIF

     !     PRINT*, ITER, OCCERROR

     IF (ABS(SHIFTCP) .GT. MAXSHIFT) THEN
        SHIFTCP = SIGN(MAXSHIFT, SHIFTCP)
     ENDIF

     CHEMPOT = CHEMPOT + SHIFTCP

#ifdef DOUBLEPREC

     IF (ITER .GE. 3 .AND. OCCERROR .LT. BREAKTOL) THEN
        BREAKLOOP = 1
     ENDIF

#elif defined(SINGLEPREC)

     IF ( ITER .GE. 3 .AND. ((OCCERROR .LT. BREAKTOL) .AND. &
          (OCCERROR .EQ. PREVERROR .OR. &
          OCCERROR .EQ. PREVERROR2 .OR. &
          OCCERROR .EQ. PREVERROR3)) .OR. &
          ITER .EQ. 25) THEN

        BREAKLOOP = 1

     ENDIF

#endif

  ENDDO

  IF (SPINON .EQ. 0) THEN
     BO = TWO*BO
  ENDIF

  RETURN

END SUBROUTINE FERMIEXPANS
