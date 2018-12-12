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

SUBROUTINE TABTEST

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE PPOTARRAY
  USE NEBLISTARRAY
  USE VIRIALARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, NEWJ, J, K, PPSEL, BREAKLOOP
  INTEGER :: PBCI, PBCJ, PBCK, COUNT
  INTEGER :: KLO, KHI
  REAL(LATTEPREC) :: RCUT2, RCUT, MAGR, MAGR2
  REAL(LATTEPREC) :: FORCE(3), DC(3), RIJ(3)
  REAL(LATTEPREC) :: GRAD, TMPE
  REAL(LATTEPREC) :: A, B, DX
  IF (EXISTERROR) RETURN

  OPEN(UNIT=47, STATUS="UNKNOWN", FILE="tabtest.dat")

  DO I = 1, 1001

     MAGR =  ONE + ONE*REAL(I-1)/1000

     KLO = 1
     KHI = PPTABLENGTH(1)

     DO WHILE (KHI - KLO .GT. 1) 

        K = (KHI + KLO)/2

        IF (PPR(K,1) .GT. MAGR) THEN
           KHI = K
        ELSE
           KLO = K
        ENDIF

     ENDDO

     DX = PPR(KHI,1) - PPR(KLO,1)

     A = (PPR(KHI, 1) - MAGR)/DX
     B = (MAGR - PPR(KLO, 1))/DX

     TMPE = A*PPVAL(KLO,1) + B*PPVAL(KHI, 1) + &
          ((A*A*A - A)*PPSPL(KLO,1) + &
          (B*B*B - B)*PPSPL(KHI,1))*(DX*DX/SIX)

     EREP = EREP + TMPE

     GRAD = (PPVAL(KHI,1) - PPVAL(KLO,1))/DX + &
          ((ONE - THREE*A*A)*PPSPL(KLO,1) + &
          (THREE*B*B - ONE)*PPSPL(KHI,1))*(DX/SIX)

     WRITE(47,*) MAGR, TMPE, GRAD

  ENDDO

  CLOSE(47)

  RETURN

END SUBROUTINE TABTEST

