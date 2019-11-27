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

SUBROUTINE DUNIVSCALE(I, J, L1, L2, MP, R, WHICHINT, SK_INTEGRAL, DSK_INTEGRAL)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE UNIVARRAY
  USE MYPRECISION

  IMPLICIT NONE
  
  INTEGER, INTENT(IN) :: I, J, L1, L2, MP
  REAL(LATTEPREC), INTENT(IN) :: R
  INTEGER :: BREAKLOOP
  INTEGER :: KLO, KHI, K, MYINTEGRAL
  REAL(LATTEPREC), INTENT(OUT) :: SK_INTEGRAL, DSK_INTEGRAL
  REAL(LATTEPREC) :: SA, SB, DX
  REAL(LATTEPREC) :: A(14),R2, R3, RMINUSR1,DPOLYNOM, RMOD
  REAL(LATTEPREC) :: POLYNOM, TMP
  REAL(LATTEPREC) :: EXP1, EXP2, EXP3, EXP4
  CHARACTER(LEN=1), INTENT(IN) :: WHICHINT
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
        
        POLYNOM = RMOD*(A(2) + RMOD*(A(3) + RMOD*(A(4) + A(5)*RMOD)))
        
        DPOLYNOM = A(2) + RMOD*(TWO*A(3) + RMOD*(THREE*A(4) + FOUR*A(5)*RMOD))
        
        SK_INTEGRAL = EXP(POLYNOM)
        
        DSK_INTEGRAL = DPOLYNOM*SK_INTEGRAL

        
     ELSEIF (R .GT. A(7) .AND. R .LT. A(8)) THEN
        
        RMINUSR1 = R - A(7)
        
        SK_INTEGRAL = A(9) + RMINUSR1*(A(10) + &
              RMINUSR1*(A(11) + RMINUSR1*(A(12) + &
              RMINUSR1*(A(13) + RMINUSR1*A(14)))))

        DSK_INTEGRAL = A(10) + RMINUSR1*(TWO*A(11) + &
             RMINUSR1*(THREE*A(12) + RMINUSR1*(FOUR*A(13) + &
             RMINUSR1*FIVE*A(14))))
        
     ELSE
        
        SK_INTEGRAL = ZERO
        DSK_INTEGRAL = ZERO
        
     END IF
     
     SK_INTEGRAL = A(1)*SK_INTEGRAL
     DSK_INTEGRAL = -A(1)*DSK_INTEGRAL

  ELSEIF (SCLTYPE .EQ. "GAUSSIAN") THEN
     
     SELECT CASE(WHICHINT)
     CASE("H") ! We're doing the H matrix build
        A = BOND(:,MYINTEGRAL)
     CASE("S") ! We're doing the S matrix build
        A = OVERL(:,MYINTEGRAL)
     END SELECT

     IF (R .LE. A(9)) THEN

        R2 = R*R
        R3 = R2*R

        EXP1 = A(1)*EXP(-A(2)*R2)
        EXP2 = A(3)*EXP(-A(4)*R2)
        EXP3 = A(5)*EXP(-A(6)*R2)
        EXP4 = A(7)*EXP(-A(8)*R2)
        
        SK_INTEGRAL = R3*EXP1 + R2*EXP2 + R*EXP3 + EXP4

        DSK_INTEGRAL = (THREE*R2 - TWO*A(2)*R3*R)*EXP1 + &
             (TWO*R - TWO*A(4)*R3)*EXP2 + &
             (ONE - TWO*A(6)*R2)*EXP3 - TWO*A(8)*R*EXP4
        
     ELSEIF (R .GT. A(9) .AND. R .LT. A(10)) THEN

!        RMINUSR1 = R - A(5)

        SK_INTEGRAL = ZERO
        DSK_INTEGRAL = ZERO

     ELSE

        SK_INTEGRAL = ZERO
        DSK_INTEGRAL = ZERO

     END IF
     
     DSK_INTEGRAL = -DSK_INTEGRAL
     
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
          
     ! Negative gradient

     IF (WHICHINT .EQ. "H") THEN
        
        SK_INTEGRAL = SA*TABH(KLO, MYINTEGRAL) + SB*TABH(KHI, MYINTEGRAL) + &
             ((SA*SA*SA - SA)*HSPL(KLO,MYINTEGRAL) + &
             (SB*SB*SB - SB)*HSPL(KHI,MYINTEGRAL))*(DX*DX/SIX)
        
        
        DSK_INTEGRAL = -((TABH(KHI,MYINTEGRAL) - TABH(KLO,MYINTEGRAL))/DX + &
             ((ONE - THREE*SA*SA)*HSPL(KLO,MYINTEGRAL) + &
             (THREE*SB*SB - ONE)*HSPL(KHI,MYINTEGRAL))*(DX/SIX))

        IF (R .GT. HCUT(MYINTEGRAL)) THEN
           SK_INTEGRAL = ZERO 
           DSK_INTEGRAL = ZERO
        ENDIF

     ELSEIF (WHICHINT .EQ. "S") THEN

        SK_INTEGRAL = SA*TABS(KLO, MYINTEGRAL) + SB*TABS(KHI, MYINTEGRAL) + &
              ((SA*SA*SA - SA)*SSPL(KLO,MYINTEGRAL) + &
              (SB*SB*SB - SB)*SSPL(KHI,MYINTEGRAL))*(DX*DX/SIX)
        
        DSK_INTEGRAL = -((TABS(KHI,MYINTEGRAL) - TABS(KLO,MYINTEGRAL))/DX + &
             ((ONE - THREE*SA*SA)*SSPL(KLO,MYINTEGRAL) + &
             (THREE*SB*SB - ONE)*SSPL(KHI,MYINTEGRAL))*(DX/SIX))
        
        IF (R .GT. SCUT(MYINTEGRAL)) THEN
           SK_INTEGRAL = ZERO
           DSK_INTEGRAL = ZERO
        ENDIF


     ENDIF
     
  ENDIF

  ! permutation symmetry
  
  IF (L1 .GT. L2 .AND. MOD(L1 + L2, 2) .NE. 0) THEN
     SK_INTEGRAL = -SK_INTEGRAL
     DSK_INTEGRAL = -DSK_INTEGRAL
  ENDIF
  
  RETURN
  
END SUBROUTINE DUNIVSCALE

