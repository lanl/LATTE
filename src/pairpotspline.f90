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

SUBROUTINE PAIRPOTSPLINE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, NEWJ, J, K, PPSEL, BREAKLOOP
  INTEGER :: PBCI, PBCJ, PBCK, COUNT
  REAL(LATTEPREC) :: RCUT2, RCUT, MAGR, MAGR2
  REAL(LATTEPREC) :: FORCE(3), DC(3), RIJ(3)
  REAL(LATTEPREC) :: GRAD, PPTMP
  REAL(LATTEPREC), EXTERNAL :: PPHEAVI

  EREP = ZERO
  FPP = ZERO
  VIRPAIR = ZERO

  DO I = 1, NATS

     ! Loop over all neighbors of I

     DO NEWJ = 1, TOTNEBPP(I)

        J = NEBPP(1, NEWJ, I)

        PBCI = NEBPP(2, NEWJ, I)
        PBCJ = NEBPP(3, NEWJ, I)
        PBCK = NEBPP(4, NEWJ, I)

        DO K = 1, NOPPS

           IF ((ATELE(I) .EQ. PPELE1(K) .AND. ATELE(J) .EQ. PPELE2(K)) &
                .OR. (ATELE(J) .EQ. PPELE1(K) .AND. &
                ATELE(I) .EQ. PPELE2(K))) THEN

              PPSEL = K
              RCUT = PPRK(1, PPSEL)
              RCUT2 = RCUT*RCUT

           ENDIF

        ENDDO

        RIJ(1) = CR(1,J) + REAL(PBCI)*BOX(1,1) + REAL(PBCJ)*BOX(2,1) + &
             REAL(PBCK)*BOX(3,1) - CR(1,I)

        RIJ(2) = CR(2,J) + REAL(PBCI)*BOX(1,2) + REAL(PBCJ)*BOX(2,2) + &
             REAL(PBCK)*BOX(3,2) - CR(2,I)

        RIJ(3) = CR(3,J) + REAL(PBCI)*BOX(1,3) + REAL(PBCJ)*BOX(2,3) + &
             REAL(PBCK)*BOX(3,3) - CR(3,I)

        MAGR2 = RIJ(1)*RIJ(1) + RIJ(2)*RIJ(2) + RIJ(3)*RIJ(3)

        IF (MAGR2 .LE. RCUT2) THEN

           MAGR = SQRT(MAGR2)

           DC = RIJ/MAGR

           GRAD = ZERO

           DO K = 1, PPNK(PPSEL)

              PPTMP = PPAK(K,PPSEL)*PPHEAVI(PPRK(K,PPSEL), MAGR)* &
                   (PPRK(K,PPSEL) - MAGR)*(PPRK(K,PPSEL) - MAGR)

              EREP = EREP + PPTMP*(PPRK(K,PPSEL) - MAGR)

              GRAD = GRAD - THREE*PPTMP

           ENDDO

           FORCE = GRAD*DC

           FPP(1,I) = FPP(1,I) + FORCE(1)
           FPP(2,I) = FPP(2,I) + FORCE(2)
           FPP(3,I) = FPP(3,I) + FORCE(3)

           VIRPAIR(1) = VIRPAIR(1) + RIJ(1)*FORCE(1)
           VIRPAIR(2) = VIRPAIR(2) + RIJ(2)*FORCE(2)
           VIRPAIR(3) = VIRPAIR(3) + RIJ(3)*FORCE(3)
           VIRPAIR(4) = VIRPAIR(4) + RIJ(1)*FORCE(2)
           VIRPAIR(5) = VIRPAIR(5) + RIJ(2)*FORCE(3)
           VIRPAIR(6) = VIRPAIR(6) + RIJ(3)*FORCE(1)

        ENDIF

     ENDDO

  ENDDO

  EREP = EREP/TWO
  VIRPAIR = VIRPAIR/TWO

  RETURN

END SUBROUTINE PAIRPOTSPLINE

FUNCTION PPHEAVI(RK,R)

  USE CONSTANTS_MOD

  IMPLICIT NONE

  REAL(LATTEPREC) :: RK, R, PPHEAVI

  IF ((RK - R) .LE. ZERO) THEN
     PPHEAVI = ZERO
  ELSE
     PPHEAVI = ONE
  ENDIF

  RETURN

END FUNCTION PPHEAVI
