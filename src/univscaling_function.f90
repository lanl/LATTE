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

! 1 = A0
! 2 = B1
! 3 = B2
! 4 = B3
! 5 = B4
! 6 = B5
! 7 = R1
! 8 = RCUT
! 9 = TAIL1
! 10 = TAIL2
! 11 = TAIL3
! 12 = TAIL4
! 13 = TAIL5
! 14 = TAIL6

FUNCTION UNIVSCALE(I, J, L1, L2, MP, R, WHICHINT)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE
  
  INTEGER :: I, J, K, L1, L2, IP1, IP2, MP, IC, MYINTEGRAL
  INTEGER :: BREAKLOOP
  INTEGER :: KLO, KHI
  REAL(LATTEPREC) :: UNIVSCALE
  REAL(LATTEPREC) :: SA, SB, DX
  REAL(LATTEPREC) :: A(14), R, RMINUSR1, RMOD, R2
  CHARACTER(LEN=1) :: WHICHINT
  CHARACTER(LEN=3) :: IGLTYPE
  IF (EXISTERROR) RETURN

  MYINTEGRAL = IGL_MAP(MP, L2, L1, J, I)

  IF (SCLTYPE .EQ. "EXP") THEN

      SELECT CASE(WHICHINT)
      CASE("H") ! We're doing the H matrix build                           
         A = BOND(:,MYINTEGRAL)
      CASE("S") ! We're doing the S matrix build                           
         A = OVERL(:,MYINTEGRAL)
      END SELECT
      
      IF (R .LE. A(7)) THEN
         
         RMOD = R - A(6)
         
         UNIVSCALE = EXP(RMOD*(A(2) + RMOD*(A(3) + RMOD*(A(4) + A(5)*RMOD))))
         
      ELSEIF (R .GT. A(7) .AND. R .LT. A(8)) THEN
         
         RMINUSR1 = R - A(7)
         
         UNIVSCALE = A(9) + RMINUSR1*(A(10) + &
              RMINUSR1*(A(11) + RMINUSR1*(A(12) + &
              RMINUSR1*(A(13) + RMINUSR1*A(14)))))
         
      ELSE
         
         UNIVSCALE = ZERO
     
      END IF
      
      UNIVSCALE = A(1)*UNIVSCALE

   ELSEIF (SCLTYPE .EQ. "GAUSSIAN") THEN
      
      SELECT CASE(WHICHINT)
      CASE("H") ! We're doing the H matrix build
         A = BOND(:,MYINTEGRAL)
      CASE("S") ! We're doing the S matrix build         
         A = OVERL(:,MYINTEGRAL)
      END SELECT

      IF (R .LE. A(9)) THEN

         R2 = R*R

         UNIVSCALE = A(1)*R2*R*EXP(-A(2)*R2) + A(3)*R2*EXP(-A(4)*R2) + &
              A(5)*R*EXP(-A(6)*R2) + A(7)*EXP(-A(8)*R2)

      ELSEIF (R .GT. A(9) .AND. R .LT. A(10)) THEN

         UNIVSCALE = ZERO

      ELSE
         
         UNIVSCALE = ZERO

      ENDIF

   ELSEIF (SCLTYPE .EQ. "TABLE") THEN

      KLO = 1
      KHI = LENTABINT(MYINTEGRAL)
      
      DO WHILE (KHI - KLO .GT. 1)
         
         K = (KHI + KLO)/2
         
         IF (TABR(K,MYINTEGRAL) .GT. R) THEN
            KHI = K
         ELSE
            KLO = K
         ENDIF
         
      ENDDO

      DX = TABR(KHI, MYINTEGRAL) - TABR(KLO,MYINTEGRAL)

      SA = (TABR(KHI, MYINTEGRAL) - R)/DX
      SB = (R - TABR(KLO, MYINTEGRAL))/DX

      IF (WHICHINT .EQ. "H") THEN 
         UNIVSCALE = SA*TABH(KLO, MYINTEGRAL) + SB*TABH(KHI, MYINTEGRAL) + &
              ((SA*SA*SA - SA)*HSPL(KLO,MYINTEGRAL) + &
              (SB*SB*SB - SB)*HSPL(KHI,MYINTEGRAL))*(DX*DX/SIX)
       
         IF (R .GT. HCUT(MYINTEGRAL)) UNIVSCALE = ZERO
  
      ELSE
         UNIVSCALE = SA*TABS(KLO, MYINTEGRAL) + SB*TABS(KHI, MYINTEGRAL) + &
              ((SA*SA*SA - SA)*SSPL(KLO,MYINTEGRAL) + &
              (SB*SB*SB - SB)*SSPL(KHI,MYINTEGRAL))*(DX*DX/SIX)
         
         IF (R .GT. SCUT(MYINTEGRAL)) UNIVSCALE = ZERO

      ENDIF

      

   ENDIF
         
  ! permutation symmetry
  
  IF (L1 .GT. L2 .AND. MOD(L1 + L2, 2) .NE. 0) UNIVSCALE = -UNIVSCALE
  
!  PRINT*, UNIVSCALE

  RETURN
  
END FUNCTION UNIVSCALE

