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

SUBROUTINE UNIVTAILCOEF(A)

  USE CONSTANTS_MOD
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I
  REAL(LATTEPREC) :: A(14)
  REAL(LATTEPREC) :: R1, R1SQ, RCUT, RMOD
  REAL(LATTEPREC) :: SCL_R1, POLYNOM, DPOLY, DDPOLY
  REAL(LATTEPREC) :: DELTA, DELTA2, DELTA3, DELTA4
  IF (EXISTERROR) RETURN

  IF ( ABS(A(1)) .LT. 0.000000000001 ) THEN

     A(9) = ZERO
     A(10) = ZERO
     A(11) = ZERO
     A(12) = ZERO
     A(13) = ZERO
     A(14) = ZERO

  ELSE

     R1 = A(7)
     RCUT = A(8)

     R1SQ = R1*R1

     RMOD = R1 - A(6)

     POLYNOM = RMOD*(A(2) + RMOD*(A(3) + RMOD*(A(4) + A(5)*RMOD)))

     SCL_R1 = EXP(POLYNOM)

     !     CALL UNIVSCALE(R1, A, SCL_R1)

     !     SCL_R1 = SCL_R1/A(1)

     DELTA = RCUT - R1

     ! Now we're using a 6th order polynomial: fitted to value, first, 
     ! and second derivatives at R1 and R_cut

     A(9) = SCL_R1

     RMOD = R1 - A(6)

     DPOLY = A(2) + TWO*A(3)*RMOD + THREE*A(4)*RMOD*RMOD + &
          FOUR*A(5)*RMOD*RMOD*RMOD

     A(10) = DPOLY*SCL_R1

     DDPOLY = TWO*A(3) + SIX*A(4)*RMOD + TWELVE*A(5)*RMOD*RMOD

     A(11) = (DPOLY*DPOLY + DDPOLY)*SCL_R1/TWO

     DELTA2 = DELTA*DELTA
     DELTA3 = DELTA2*DELTA
     DELTA4 = DELTA3*DELTA

     A(12) = (MINUSONE/DELTA3)*(THREE*A(11)*DELTA2 + &
          SIX*A(10)*DELTA + TEN*A(9))

     A(13) = (ONE/DELTA4)*(THREE*A(11)*DELTA2 + &
          EIGHT*A(10)*DELTA + FIFTEEN*A(9))

     A(14) = (MINUSONE/(TEN*DELTA3)) * &
          (SIX*A(13)*DELTA2 + THREE*A(12)*DELTA + A(11))

  ENDIF

  RETURN

END SUBROUTINE UNIVTAILCOEF
