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

SUBROUTINE SP2DBCSR

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


  INTEGER :: I, J, ITER
  INTEGER :: BREAKLOOP
  INTEGER :: PUR
  REAL(LATTEPREC) :: TRX, OCC, TRX2, LIMDIFF
  REAL(LATTEPREC) :: GERSHFACT
  REAL(LATTEPREC) :: IDEMPERR, TRXOLD
  REAL(LATTEPREC) :: IDEMPERR1, IDEMPERR2
  REAL(LATTEPREC), PARAMETER :: IDEMTOL = 1.0D-14


  IDEMPERR = ZERO
  IDEMPERR1 = ZERO
  IDEMPERR2 = ZERO

  !
  ! We're also using Niklasson's scheme to determine convergence
  !

  OCC = BNDFIL*FLOAT(HDIM)

  ! Using intrinsics is probably better than coding this ourselves

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

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, MINUSONE, &
          BO, HDIM, BO, HDIM, ONE, X2, HDIM)


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
10   FORMAT(I4, 2G30.18)
     IF (SP2CONV .EQ. "REL" .AND. ITER .GE. MINSP2ITER &
          .AND. (IDEMPERR2 .LE. IDEMPERR .OR. &
          IDEMPERR .LT. IDEMTOL)) BREAKLOOP = 1

     !	IF (ITER .EQ. 30) BREAKLOOP=1
     IF (SP2CONV .EQ. "ABS" .AND. ABS(LIMDIFF) .LT. IDEMTOL) BREAKLOOP = 1

  ENDDO

  IF (ITER .EQ. 100) THEN
     CALL PANIC
     CALL ERRORS("sp2dbcsr","SP2 purification is not converging: STOP!")
  ENDIF

  BO = TWO*BO


  RETURN

END SUBROUTINE SP2DBCSR
