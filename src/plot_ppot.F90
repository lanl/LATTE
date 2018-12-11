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

SUBROUTINE PLOTPPOT

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, K, II, PPSEL
  REAL(LATTEPREC) :: JR1, JRCUT, R1, RCUT2, RCUT
  REAL(LATTEPREC) :: PHI, FORCE
  REAL(LATTEPREC) :: MYR, MYR2, MYR3, MYR4, MAGR2, MAGR, R6
  REAL(LATTEPREC), ALLOCATABLE :: MYPAIRPOTS(:,:), MYF(:,:)
  REAL(LATTEPREC) :: DC(3), DPHI(3), EXPTMP
  REAL(LATTEPREC) :: POLYNOM, DPOLYNOM
  CHARACTER(LEN=50) :: FLNM
  IF (EXISTERROR) RETURN

  ! There are 4 different parts to the pair potential:
  !
  ! 1) Short range repulsion fitting to give bond lengths etc
  ! 2) The joining function from JOINR1 TO JOINRCUT
  ! 3) The vdW-type pair potential from JOINCUT to PPR1
  ! 4) The final cut off tail from PPR1 TO PPRCUT
  !

  OPEN(UNIT=56, STATUS="UNKNOWN", FILE = "vdW_scaling.dat")
  OPEN(UNIT=57, STATUS="UNKNOWN", FILE = "vdW_F_scaling.dat")

  ALLOCATE(MYPAIRPOTS(1000,NOPPS), MYF(1000,NOPPS))

  MYPAIRPOTS = ZERO

  DC(1) = ONE
  DC(2) = ZERO
  DC(3) = ZERO

  DO PPSEL = 1, NOPPS

     R1 = POTCOEF(9, PPSEL)
     RCUT = POTCOEF(10,PPSEL)

     DO II = 1, 1000

        MAGR = 0.7D0 + 4.3D0*REAL(II-1)/THOUSAND

        !        MAGR = MAGR - POTCOEF(6,PPSEL)

        MAGR2 = MAGR*MAGR

        IF (MAGR .LT. R1) THEN

           MAGR = MAGR - POTCOEF(6,PPSEL)

           !           CALL DUNIVSCALE(MAGR, POTCOEF(:,PPSEL), DC, PHI, DPHI)

           POLYNOM = MAGR*(POTCOEF(2,PPSEL) + MAGR*(POTCOEF(3,PPSEL) + &
                MAGR*(POTCOEF(4,PPSEL) + MAGR*POTCOEF(5,PPSEL))))

           PHI = POTCOEF(1,PPSEL)*EXP(POLYNOM)

           DPOLYNOM = POTCOEF(2,PPSEL) + MAGR*(TWO*POTCOEF(3,PPSEL) + &
                MAGR*(THREE*POTCOEF(4,PPSEL) + &
                FOUR*POTCOEF(5,PPSEL)*MAGR))

           DPHI = -DC*PHI*DPOLYNOM           

           !           EXPTMP = POTCOEF(6,PPSEL)*EXP(POTCOEF(7,PPSEL) * &
           !                (MAGR - POTCOEF(8,PPSEL)))

           EXPTMP = ZERO

           MYPAIRPOTS(II, PPSEL) = PHI + EXPTMP

           !           MYPAIRPOTS(II, PPSEL) =  PHI + EXPTMP - &
           !                POTCOEF(8,PPSEL)/(MAGR2*MAGR2*MAGR2)

           !           MYF(II,PPSEL) = -DPHI(1) + DC(1)*(POTCOEF(7,PPSEL)*EXPTMP + &
           !                SIX*POTCOEF(8,PPSEL)/(MAGR*MAGR2*MAGR2*MAGR2))

        ELSEIF (MAGR .GE. R1 .AND. MAGR .LT. RCUT) THEN

           MYR = MAGR - R1

           MYPAIRPOTS(II, PPSEL) = POTCOEF(11,PPSEL) + &
                MYR*(POTCOEF(12,PPSEL) + MYR*(POTCOEF(13,PPSEL) + &
                MYR*(POTCOEF(14,PPSEL) + MYR*(POTCOEF(15,PPSEL) + &
                MYR*POTCOEF(16,PPSEL)))))

           MYF(II,PPSEL) = POTCOEF(12,PPSEL) +  MYR*(TWO*POTCOEF(13,PPSEL) + &
                MYR*(THREE*POTCOEF(14,PPSEL) + &
                MYR*(FOUR*POTCOEF(15,PPSEL) + &
                MYR*FIVE*POTCOEF(16,PPSEL))))

        ENDIF

     ENDDO

  ENDDO

  DO II = 1, 1000

     MAGR = 0.7D0 + 4.3D0*REAL(II-1)/THOUSAND

     WRITE(56, 10) MAGR, &
          (MYPAIRPOTS(II,I), I = 1, NOPPS)

     WRITE(57, 10) MAGR , &
          (MYF(II,I), I = 1, NOPPS)
  ENDDO

10 FORMAT(100G18.9)

  CLOSE(56)
  CLOSE(57)

  DEALLOCATE(MYPAIRPOTS, MYF)

  RETURN

END SUBROUTINE PLOTPPOT

