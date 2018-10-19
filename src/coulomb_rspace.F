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

SUBROUTINE COULOMBRSPACE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE COULOMBARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, NEWJ
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: RIJ(3), MAGR, MAGR2, DC(3)
  REAL(LATTEPREC) :: CA, TI, TJ, FORCE
  REAL(LATTEPREC) :: TI2, TI3, TI4, TI6, TJ2, TJ3, TJ4, TJ6
  REAL(LATTEPREC) :: EXPTI, EXPTJ
  REAL(LATTEPREC) :: TI2MTJ2, TJ2MTI2
  REAL(LATTEPREC) :: SA, SB, SC, SD, SE, SF
  REAL(LATTEPREC) :: SSA, SSB, SSC, SSD, SSE
  REAL(LATTEPREC), PARAMETER :: C1 = -1.26551223, C2 = 1.00002368
  REAL(LATTEPREC), PARAMETER :: C3 = 0.37409196, C4 = 0.09678418
  REAL(LATTEPREC), PARAMETER :: C5 = -0.18628806, C6 = 0.27886807
  REAL(LATTEPREC), PARAMETER :: C7 = -1.13520398, C8 = 1.48851587
  REAL(LATTEPREC), PARAMETER :: C9 = -0.82215223, C10 = 0.17087277 
  REAL(LATTEPREC) :: Z, T, NUMREP_ERFC
  IF (EXISTERROR) RETURN
  !  REAL(LATTEPREC), EXTERNAL :: NUMREP_ERFC

  FCOUL = ZERO
  COULOMBV = ZERO

  VIRCOUL = ZERO

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(NATS, BOX, CR, COULCUT2, TOTNEBCOUL, NEBCOUL, HUBBARDU, ELEMPOINTER, DELTAQ, COULOMBV ) &
  !$OMP SHARED(CALPHA, CALPHA2, SQRTPI, FCOUL, KECONST, TFACT) &
  !$OMP PRIVATE(I, J, NEWJ, PBCI, PBCJ, PBCK, RIJ, MAGR2, MAGR,TI,TI2,TI3,TI4,TI6,SSA ) &
  !$OMP PRIVATE(SSB, SSC, SSD, SSE, SA, SB, SC, SD, SE, SF ) &
  !$OMP PRIVATE(TI2MTJ2, TJ2MTI2, EXPTI, EXPTJ, TJ, TJ2, TJ3, TJ4, TJ6, FORCE, CA, DC, Z, T, NUMREP_ERFC) &
  !$OMP REDUCTION(+:VIRCOUL)  
  DO I = 1, NATS

     TI = TFACT*HUBBARDU(ELEMPOINTER(I))

     TI2 = TI*TI
     TI3 = TI2*TI
     TI4 = TI2*TI2
     TI6 = TI4*TI2

     SSA = TI 
     SSB = TI3/FORTYEIGHT
     SSC = THREE*TI2/SIXTEEN
     SSD = ELEVEN*TI/SIXTEEN
     SSE = ONE

     DO NEWJ = 1, TOTNEBCOUL(I)

        J = NEBCOUL(1, NEWJ, I)
        PBCI = NEBCOUL(2, NEWJ, I)
        PBCJ = NEBCOUL(3, NEWJ, I)
        PBCK = NEBCOUL(4, NEWJ, I)

        RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1) - CR(1,I)

        RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2) - CR(2,I)

        RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3) - CR(3,I)

        MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

        IF (MAGR2 .LE. COULCUT2) THEN

           TJ = TFACT*HUBBARDU(ELEMPOINTER(J)) 

           MAGR = SQRT(MAGR2)

           !
           ! Direction cosines
           !

           DC = RIJ/MAGR

           ! Using Numerical Recipes ERFC

           Z = ABS(CALPHA*MAGR)

           T = ONE/(ONE + HALF*Z)

           NUMREP_ERFC = T * EXP(-Z*Z + C1 + T*(C2 + T*(C3 + T*(C4 + T*(C5 + &
                T*(C6 + T*(C7 + T*(C8 + T*(C9 + T*C10)))))))))

           IF ( CALPHA*MAGR .LT. ZERO ) NUMREP_ERFC = TWO - NUMREP_ERFC

           CA = NUMREP_ERFC/MAGR

           COULOMBV(I) = COULOMBV(I) + DELTAQ(J)*CA

           CA = CA + TWO*CALPHA*EXP( -CALPHA2*MAGR2 )/SQRTPI

           FORCE = -KECONST*DELTAQ(I)*DELTAQ(J)*CA/MAGR

           EXPTI = EXP( -TI*MAGR )

           IF (ELEMPOINTER(I) .EQ. ELEMPOINTER(J)) THEN

              COULOMBV(I) = COULOMBV(I) - DELTAQ(J)*EXPTI* &
                   (SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR)

              FORCE = FORCE + (KECONST*DELTAQ(I)*DELTAQ(J)*EXPTI) * &
                   ((SSE/MAGR2 - TWO*SSB*MAGR - SSC) + &
                   SSA*(SSB*MAGR2 + SSC*MAGR + SSD + SSE/MAGR))

           ELSE

              TJ2 = TJ*TJ
              TJ3 = TJ2*TJ
              TJ4 = TJ2*TJ2
              TJ6 = TJ4*TJ2

              EXPTJ = EXP( -TJ*MAGR )

              TI2MTJ2 = TI2 - TJ2
              TJ2MTI2 = -TI2MTJ2

              SA = TI
              SB = TJ4*TI/(TWO * TI2MTJ2 * TI2MTJ2)
              SC = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2 * TI2MTJ2 * TI2MTJ2)

              SD = TJ
              SE = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
              SF = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2 * TJ2MTI2 * TJ2MTI2)

              COULOMBV(I) = COULOMBV(I) - (DELTAQ(J) * &
                   (EXPTI*(SB - (SC/MAGR)) + EXPTJ*(SE - (SF/MAGR))))

              FORCE = FORCE + KECONST*DELTAQ(I)*DELTAQ(J) * &
                   ((EXPTI * (SA*(SB - (SC/MAGR)) - (SC/MAGR2))) + &
                   (EXPTJ * (SD*(SE - (SF/MAGR)) - (SF/MAGR2))))

           ENDIF

           FCOUL(1,I) = FCOUL(1,I) + DC(1)*FORCE
           FCOUL(2,I) = FCOUL(2,I) + DC(2)*FORCE
           FCOUL(3,I) = FCOUL(3,I) + DC(3)*FORCE

           VIRCOUL(1) = VIRCOUL(1) + RIJ(1)*DC(1)*FORCE
           VIRCOUL(2) = VIRCOUL(2) + RIJ(2)*DC(2)*FORCE
           VIRCOUL(3) = VIRCOUL(3) + RIJ(3)*DC(3)*FORCE
           VIRCOUL(4) = VIRCOUL(4) + RIJ(1)*DC(2)*FORCE
           VIRCOUL(5) = VIRCOUL(5) + RIJ(2)*DC(3)*FORCE
           VIRCOUL(6) = VIRCOUL(6) + RIJ(3)*DC(1)*FORCE

        ENDIF

     ENDDO

  ENDDO
  !$OMP END PARALLEL DO


  COULOMBV = KECONST * COULOMBV
  VIRCOUL = HALF * VIRCOUL

  RETURN

END SUBROUTINE COULOMBRSPACE


! Taken from: Numerical Recipes in Fortran 77. The Art of Scientific
! Computing (ISBN 0-521-43064-X).

!FUNCTION NUMREP_ERFC( X )

!  USE MYPRECISION

! IMPLICIT NONE

!  REAL(LATTEPREC) :: NUMREP_ERFC, X
!  REAL(LATTEPREC) :: T, Z
!  REAL(LATTEPREC), PARAMETER :: C1 = -1.26551223, C2 = 1.00002368
! REAL(LATTEPREC), PARAMETER :: C3 = 0.37409196, C4 = 0.09678418
!  REAL(LATTEPREC), PARAMETER :: C5 = -0.18628806, C6 = 0.27886807
!  REAL(LATTEPREC), PARAMETER :: C7 = -1.13520398, C8 = 1.48851587
!  REAL(LATTEPREC), PARAMETER :: C9 = -0.82215223, C10 = 0.17087277 

!  Z = ABS(X)

!  T = ONE/(ONE + HALF*Z)

!  NUMREP_ERFC = T * EXP(-Z*Z + C1 + T*(C2 + T*(C3 + T*(C4 + T*(C5 + &
!       T*(C6 + T*(C7 + T*(C8 + T*(C9 + T*C10)))))))))

!  IF ( X .LT. ZERO ) NUMREP_ERFC = TWO - NUMREP_ERFC

!  RETURN

!END FUNCTION NUMREP_ERFC


