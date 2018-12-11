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

SUBROUTINE PAIRPOTNONEB

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, NEWJ, J, K, PPSEL, BREAKLOOP
  INTEGER :: PBCI, PBCJ, PBCK
  REAL(LATTEPREC) :: JR1, JRCUT, R1, RCUT2, RCUT
  REAL(LATTEPREC) :: FORCE, DC(3), RIJ(3)
  REAL(LATTEPREC) :: MYR, MYR2, MYR3, MYR4, MAGR2, MAGR
  REAL(LATTEPREC) :: UNIVPHI, JOINPHI, VDWPHI, CUTPHI, TMP
  REAL(LATTEPREC) :: VIRUNIV(6), VIRJOIN(6), VIRVDW(6), VIRCUT(6)
  REAL(LATTEPREC) :: FUNIV(3), FJOIN(3), FVDW(3), FCUT(3)
  REAL(LATTEPREC) :: PHI, DPHI(3), EXPTMP, R6, FTMP(3)
  REAL(LATTEPREC) :: POLYNOM, DPOLYNOM
  IF (EXISTERROR) RETURN

  ! 
  ! In this subroutine we add contributions in a strange way to ensure
  ! numerical accuracy when switching between single and double precision.
  ! If we don't do this, we get errors associated with adding very small
  ! numbers to very large ones, and energies can be off by 0.01% or more.
  !

  !
  ! There are 4 different parts to the pair potential:
  !
  ! 1) Short range repulsion fitting to give bond lengths etc
  ! 2) The joining function from JOINR1 TO JOINRCUT
  ! 3) The vdW-type pair potential from JOINCUT to PPR1
  ! 4) The final cut off tail from PPR1 TO PPRCUT
  !

  UNIVPHI = ZERO
  CUTPHI = ZERO

  VIRUNIV = ZERO
  VIRCUT = ZERO

  IF (PPOTON .EQ. 1) THEN

     DO I = 1, NATS

        FUNIV = ZERO
        FCUT = ZERO

        ! Loop over all neighbors of I

        DO J = 1, NATS

           !        DO NEWJ = 1, TOTNEBPP(I)

           !           J = NEBPP(1, NEWJ, I)

           !           PBCI = NEBPP(2, NEWJ, I)
           !           PBCJ = NEBPP(3, NEWJ, I)
           !           PBCK = NEBPP(4, NEWJ, I)

           DO K = 1, NOPPS

              IF ((ATELE(I) .EQ. PPELE1(K) .AND. ATELE(J) .EQ. PPELE2(K)) &
                   .OR. (ATELE(J) .EQ. PPELE1(K) .AND. &
                   ATELE(I) .EQ. PPELE2(K))) THEN

                 PPSEL = K

                 R1 = POTCOEF(9,PPSEL)
                 RCUT = POTCOEF(10,PPSEL)
                 RCUT2 = RCUT*RCUT

              ENDIF

           ENDDO

           RIJ(1) = CR(1,J) - CR(1,I)
           RIJ(2) = CR(2,J) - CR(2,I)
           RIJ(3) = CR(3,J) - CR(3,I)

           !           RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
           !                REAL(PBCK)*BOX(3,1) - CR(1,I)

           !           RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
           !                REAL(PBCK)*BOX(3,2) - CR(2,I)

           !           RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
           !                REAL(PBCK)*BOX(3,3) - CR(3,I)

           MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

           IF (MAGR2 .LE. RCUT2 .AND. J .NE. I) THEN

              MAGR = SQRT(MAGR2)

              ! Direction cosines

              DC = RIJ/MAGR

              IF (MAGR .LT. R1) THEN

                 !                 CALL DUNIVSCALE(MAGR, POTCOEF(:,PPSEL), DC, PHI, DPHI)

                 POLYNOM = MAGR*(POTCOEF(2,PPSEL) + MAGR*(POTCOEF(3,PPSEL) + &
                      MAGR*(POTCOEF(4,PPSEL) + MAGR*POTCOEF(5,PPSEL))))

                 PHI = POTCOEF(1,PPSEL)*EXP(POLYNOM)

                 DPOLYNOM = POTCOEF(2,PPSEL) + MAGR*(TWO*POTCOEF(3,PPSEL) + &
                      MAGR*(THREE*POTCOEF(4,PPSEL) + &
                      FOUR*POTCOEF(5,PPSEL)*MAGR))

                 DPHI = -DC*PHI*DPOLYNOM

                 ! Hack!
                 EXPTMP = POTCOEF(6,PPSEL)*&
                      EXP( POTCOEF(7,PPSEL) * (MAGR - POTCOEF(8,PPSEL)) )
                 !                 R6 = MAGR2*MAGR2*MAGR2

                 !                 UNIVPHI = UNIVPHI + PHI + EXPTMP - POTCOEF(8,PPSEL)/R6

                 UNIVPHI = UNIVPHI + PHI + EXPTMP

                 FTMP = DC*POTCOEF(7,PPSEL)*EXPTMP
                 !                 FTMP = DC*(POTCOEF(7,PPSEL)*EXPTMP + &
                 !                      SIX*POTCOEF(8,PPSEL)/(MAGR*R6))

                 FUNIV = FUNIV - DPHI + FTMP

                 VIRUNIV(1) = VIRUNIV(1) - RIJ(1)*(DPHI(1) - FTMP(1))
                 VIRUNIV(2) = VIRUNIV(2) - RIJ(2)*(DPHI(2) - FTMP(2))
                 VIRUNIV(3) = VIRUNIV(3) - RIJ(3)*(DPHI(3) - FTMP(3))
                 VIRUNIV(4) = VIRUNIV(4) - RIJ(1)*(DPHI(2) - FTMP(2))
                 VIRUNIV(5) = VIRUNIV(5) - RIJ(2)*(DPHI(3) - FTMP(3))
                 VIRUNIV(6) = VIRUNIV(6) - RIJ(3)*(DPHI(1) - FTMP(1))

              ELSE

                 MYR = MAGR - R1

                 CUTPHI =  CUTPHI + POTCOEF(11,PPSEL) + &
                      MYR*(POTCOEF(12,PPSEL) + MYR*(POTCOEF(13,PPSEL) + &
                      MYR*(POTCOEF(14,PPSEL) + MYR*(POTCOEF(15,PPSEL) + &
                      MYR*POTCOEF(16,PPSEL)))))

                 FORCE = POTCOEF(12,PPSEL)  + MYR*(TWO*POTCOEF(13,PPSEL) + &
                      MYR*(THREE*POTCOEF(14,PPSEL) + &
                      MYR*(FOUR*POTCOEF(15,PPSEL) + &
                      MYR*FIVE*POTCOEF(16,PPSEL))))

                 FCUT = FCUT + DC*FORCE

                 VIRCUT(1) = VIRCUT(1) + RIJ(1)*DC(1)*FORCE
                 VIRCUT(2) = VIRCUT(2) + RIJ(2)*DC(2)*FORCE
                 VIRCUT(3) = VIRCUT(3) + RIJ(3)*DC(3)*FORCE
                 VIRCUT(4) = VIRCUT(4) + RIJ(1)*DC(2)*FORCE
                 VIRCUT(5) = VIRCUT(5) + RIJ(2)*DC(3)*FORCE
                 VIRCUT(6) = VIRCUT(6) + RIJ(3)*DC(1)*FORCE 

              ENDIF

           ENDIF

        ENDDO

        FPP(1,I) = FUNIV(1) + FCUT(1)
        FPP(2,I) = FUNIV(2) + FCUT(2)
        FPP(3,I) = FUNIV(3) + FCUT(3)

     ENDDO

     EREP = HALF*(UNIVPHI + CUTPHI)

     !  PRINT*, "EREP ", EREP
     VIRPAIR = HALF*(VIRUNIV + VIRCUT)

  ELSE

     FPP = ZERO
     EREP = ZERO
     VIRPAIR = ZERO

  ENDIF

  RETURN

END SUBROUTINE PAIRPOTNONEB

