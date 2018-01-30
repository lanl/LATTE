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

SUBROUTINE SP2T

  !
  ! This subroutine implements the finite temperature version of
  ! Niklasson's SP2 purification scheme. It will give us something that
  ! looks very similar to a Fermi-Dirac distribution (although it's 
  ! not exactly the FD distribution).
  !

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PUREARRAY
  USE SPINARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER, II
  REAL(LATTEPREC) :: TRX2, TRX
  REAL(LATTEPREC) :: OCC, LAMBDA
  IF (EXISTERROR) RETURN

  OCC = BNDFIL*FLOAT(HDIM)

  IF (SPINON .EQ. 0) THEN

     !
     ! Rescale H such that eigenvalues are in the range [0:1] based
     ! on Gershgorin circles
     !

     DO I = 1, HDIM
        DO J = I, HDIM

           IF (I .EQ. J) THEN

              BO(J,I) = (MAXEVAL - H(J,I))/MAXMINUSMIN

           ELSEIF (I .NE. J) THEN

              BO(J,I) = (ZERO - H(J,I))/MAXMINUSMIN
              BO(I,J) = BO(J,I)

           ENDIF

        ENDDO
     ENDDO

     ITER = 0

     DO II = 1, NORECS

        ITER = ITER + 1

        !
        ! Calculating X*X
        !

#ifdef DOUBLEPREV
        CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
             BO, HDIM, BO, HDIM, 0.0D0, X2, HDIM)
#elif defined(SINGLEPREC)
        CALL SGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0, &
             BO, HDIM, BO, HDIM, 0.0, X2, HDIM)
#endif

        TRX2 = ZERO
        TRX = ZERO

        DO I = 1, HDIM

           TRX2 = TRX2 + X2(I,I)
           TRX = TRX + BO(I,I)

        ENDDO

        IF (ABS(TRX2 - OCC) .LT. ABS(TWO*TRX - TRX2 - OCC)) THEN

           BO = X2

        ELSE

           BO = TWO*BO - X2

        ENDIF

     ENDDO

     !
     ! Now apply a shift such that we get the correct occupation after
     ! only NOREC applications of the purification algorithm
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

     LAMBDA = (OCC - TRX)/(TRX - TRX2)

     !
     ! Here we also double the density matrix to get the 
     ! bond-order in non-spin polarized calculations
     !

     BO = TWO*( BO + LAMBDA*( BO - X2 ))

  ELSE

     !
     ! Remapping eigenvalues of H
     !

     DO I = 1, HDIM
        DO J = I, HDIM

           IF (I .EQ. J) THEN

              RHOUP(J,I) = (MAXEVAL - HUP(J,I))/MAXMINUSMIN

              RHODOWN(J,I) = (MAXEVAL - HDOWN(J,I)) / &
                   MAXMINUSMIN

           ELSEIF (I .NE. J) THEN

              RHOUP(J,I) = (ZERO - HUP(J,I))/MAXMINUSMIN
              RHOUP(I,J) = RHOUP(J,I)

              RHODOWN(J,I) = (ZERO - HDOWN(J,I))/MAXMINUSMIN
              RHODOWN(I,J) = RHODOWN(J,I)

           ENDIF

        ENDDO
     ENDDO

     ITER = 0

     DO II = 1, NORECS

        ITER = ITER + 1

        !
        ! Calculating Xup*Xup and Xdown*Xdown
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

           RHOUP = X2UP
           RHODOWN = X2DOWN

        ELSE

           RHOUP = TWO*RHOUP - X2UP
           RHODOWN = TWO*RHODOWN - X2DOWN

        ENDIF

     ENDDO

     !
     ! Apply a shift such that we get the correct trace after NORECS
     ! SP2 recursions
     !

#ifdef DOUBLEPREC

     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHOUP, HDIM, RHOUP, HDIM, 0.0D0, X2UP, HDIM)
     CALL DGEMM('N', 'N', HDIM, HDIM, HDIM, 1.0D0, &
          RHODOwN, HDIM, RHODOWN, HDIM, 0.0D0, X2DOWN, HDIM)

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

     IF (ABS(TOTNE - TRX) .GT. 1.0D-8) THEN

        LAMBDA = (TOTNE - TRX)/(TRX - TRX2)

        !     print*, "TOTNE - TRX = ", TOTNE - TRX
        !     print*, "TRX - TRX2 = ", TRX - TRX2

        RHOUP = RHOUP + LAMBDA*(RHOUP - X2UP)
        RHODOWN = RHODOWN + LAMBDA*(RHODOWN - X2DOWN)   

     ENDIF

  ENDIF

  RETURN

END SUBROUTINE SP2T
