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

  INTEGER :: I, J, K, NEWJ, JJ, II
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: RIJ(3), MAGR, MAGR2, DC(3)
  REAL(LATTEPREC) :: CA, TI, TJ, FORCE
  REAL(LATTEPREC) :: TI2, TI3, TI4, TI6, TJ2, TJ3, TJ4, TJ6
  REAL(LATTEPREC) :: EXPTI, EXPTJ
  REAL(LATTEPREC) :: TI2MTJ2, TJ2MTI2
  REAL(LATTEPREC) :: Z, T, NUMREP_ERFC
  REAL(LATTEPREC),SAVE, ALLOCATABLE :: TS(:), SSA(:), SSB(:), SSC(:)
  REAL(LATTEPREC),SAVE, ALLOCATABLE :: SSD(:), SSE(:)
  REAL(LATTEPREC),SAVE, ALLOCATABLE :: SB(:,:), SC(:,:)
  REAL(LATTEPREC),SAVE, ALLOCATABLE :: SE(:,:), SF(:,:)

  IF (EXISTERROR) RETURN

  FCOUL = ZERO
  COULOMBV = ZERO
  VIRCOUL = ZERO

  IF (.NOT. ALLOCATED(TS)) THEN
     ALLOCATE(TS(SIZE(HUBBARDU,1)))
     ALLOCATE(SSA(SIZE(HUBBARDU,1)))
     ALLOCATE(SSB(SIZE(HUBBARDU,1)))
     ALLOCATE(SSC(SIZE(HUBBARDU,1)))
     ALLOCATE(SSD(SIZE(HUBBARDU,1)))
     ALLOCATE(SSE(SIZE(HUBBARDU,1)))

     ALLOCATE(SB(SIZE(HUBBARDU,1),SIZE(HUBBARDU,1)))
     ALLOCATE(SC(SIZE(HUBBARDU,1),SIZE(HUBBARDU,1)))
     ALLOCATE(SE(SIZE(HUBBARDU,1),SIZE(HUBBARDU,1)))
     ALLOCATE(SF(SIZE(HUBBARDU,1),SIZE(HUBBARDU,1)))

     DO I = 1, SIZE(HUBBARDU,1)
        TS(I) = TFACT*HUBBARDU(I)
        TI2 = TS(I)**2
        TI3 = TI2*TS(I)
        TI4 = TI2*TI2
        TI6 = TI4*TI2
        SSA(I) = TS(I)
        SSB(I) = TI3/FORTYEIGHT
        SSC(I) = THREE*TI2/SIXTEEN
        SSD(I) = ELEVEN*TS(I)/SIXTEEN
        SSE(I) = ONE
        DO J = 1, SIZE(HUBBARDU,1)
           TJ = TFACT*HUBBARDU(J)
           TJ2 = TJ*TJ
           TJ3 = TJ2*TJ
           TJ4 = TJ2*TJ2
           TJ6 = TJ4*TJ2
           TI2MTJ2 = TI2 - TJ2
           TJ2MTI2 = -TI2MTJ2
           SB(I,J) = TJ4*TS(I)/(TWO * TI2MTJ2 * TI2MTJ2)
           SC(I,J) = (TJ6 - THREE*TJ4*TI2)/(TI2MTJ2 * TI2MTJ2 * TI2MTJ2)
           SE(I,J) = TI4*TJ/(TWO * TJ2MTI2 * TJ2MTI2)
           SF(I,J) = (TI6 - THREE*TI4*TJ2)/(TJ2MTI2 * TJ2MTI2 * TJ2MTI2)
        ENDDO
     ENDDO

  END IF

  !$OMP PARALLEL DO DEFAULT(NONE) &
  !$OMP SHARED(NATS, BOX, CR, COULCUT2, TOTNEBCOUL, NEBCOUL, HUBBARDU, ELEMPOINTER, DELTAQ, COULOMBV ) &
  !$OMP SHARED(CALPHA, CALPHA2, SQRTPI, FCOUL, KECONST, TFACT, TS) &
  !$OMP SHARED(SSA, SSB, SSC, SSD, SSE, SB, SC, SE, SF) &
  !$OMP PRIVATE(I, J, NEWJ, PBCI, PBCJ, PBCK, RIJ, MAGR2, MAGR, TI, II, JJ) &
  !$OMP PRIVATE(EXPTI, EXPTJ, TJ, FORCE, CA, DC, Z) &
  !$OMP REDUCTION(+:VIRCOUL)
  DO I = 1, NATS

     II = ELEMPOINTER(I)
     TI = TS(II)

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

           JJ = ELEMPOINTER(J)
           TJ = TFACT*HUBBARDU(JJ)

           MAGR = SQRT(MAGR2)

           !
           ! Direction cosines
           !

           DC = RIJ/MAGR

           Z = ABS(CALPHA*MAGR)

           CA = ERFC(Z)/MAGR

           COULOMBV(I) = COULOMBV(I) + DELTAQ(J)*CA

           CA = CA + TWO*CALPHA*EXP( -CALPHA2*MAGR2 )/SQRTPI

           FORCE = -KECONST*DELTAQ(I)*DELTAQ(J)*CA/MAGR

           EXPTI = EXP( -TI*MAGR )

           IF (ELEMPOINTER(I) .EQ. ELEMPOINTER(J)) THEN

              COULOMBV(I) = COULOMBV(I) - DELTAQ(J)*EXPTI* &
                   (SSB(JJ)*MAGR2 + SSC(JJ)*MAGR &
                   & + SSD(JJ) + SSE(JJ)/MAGR)

              FORCE = FORCE + (KECONST*DELTAQ(I)*DELTAQ(J)*EXPTI) * &
                   ((SSE(JJ)/MAGR2 - TWO*SSB(JJ)*MAGR &
                   & - SSC(JJ)) + SSA(JJ)*(SSB(JJ)*MAGR2 + &
                   SSC(JJ)*MAGR + SSD(JJ) + SSE(JJ)/MAGR))

           ELSE

              TJ = TS(JJ)
              EXPTJ = EXP( -TJ*MAGR )

              COULOMBV(I) = COULOMBV(I) - (DELTAQ(J) * &
                   & (EXPTI*(SB(II,JJ) - (SC(II,JJ)/MAGR)) + &
                   & EXPTJ*(SE(II,JJ) - (SF(II,JJ)/MAGR))))

              FORCE = FORCE + KECONST*DELTAQ(I)*DELTAQ(J) * &
                   & ((EXPTI * (TI*(SB(II,JJ) - (SC(II,JJ)/MAGR)) - &
                   & (SC(II,JJ)/MAGR2))) + (EXPTJ * (TJ*(SE(II,JJ) - &
                   & (SF(II,JJ)/MAGR)) - (SF(II,JJ)/MAGR2))))

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
