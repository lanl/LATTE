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

SUBROUTINE VDWTAILCOEF

  USE CONSTANTS_MOD
  USE PPOTARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, PPID
  REAL(LATTEPREC) :: A(17)
  REAL(LATTEPREC) :: DELTA, DELTA2, DELTA3, DELTA4
  REAL(LATTEPREC) :: X, Y, Z, MYR, R6
  REAL(LATTEPREC) :: R1, R1SQ, RCUT
  REAL(LATTEPREC) :: SCL_R1, POLY, DPOLY, DDPOLY, EXPTMP
  IF (EXISTERROR) RETURN

  ! POTCOEF:
  ! 1  2  3  4  5  6  7  8 9  10   11 12 13 14 15 16
  ! A0 A1 A2 A3 A4 A5 A6 C R1 RCUT B1 B2 B3 B4 B5 B6

  ! PHI = A0*EXP(A1*X + A2*X^2 + A3*X^3 + A4*X^4) + A5*EXP(A6*X) - C/X^6

  ! Now lets' try this:

  ! Phi A0*EXP(A1*(X - A5) + A2*(X - A5)^2 + A3*(X-A5)^3 + A4*(X-A5)^4)



  !
  ! The cut-offs and joining functions look like this:
  !
  ! t(R) = B1 + B2*(R - R1) + B3*(R - R1)^2 + B4*(R - R1)^3 + 
  !          B5*(R - R1)^4 + B6*(R - R1)^5
  !

  DO PPID = 1, NOPPS

     ! Join first : need values for pair potential at R1

     R1 = POTCOEF(9,PPID) - POTCOEF(6,PPID)
     RCUT = POTCOEF(10,PPID)
     R1SQ = R1*R1
     !     R6 = R1SQ * R1SQ * R1SQ

     !     CALL UNIVSCALE(R1, POTCOEF(:,PPID), SCL_R1)

     POLY = R1*(POTCOEF(2,PPID) + R1*(POTCOEF(3,PPID) + &
          R1*(POTCOEF(4,PPID) + R1*POTCOEF(5,PPID))))

     SCL_R1 = POTCOEF(1,PPID)*EXP(POLY)     

     !     EXPTMP = POTCOEF(6,PPID)*EXP(POTCOEF(7,PPID)*(R1 - POTCOEF(8,PPID)))

     EXPTMP = ZERO

     POTCOEF(11,PPID) = SCL_R1 + EXPTMP ! - POTCOEF(8,PPID)/R6

     DPOLY = POTCOEF(2,PPID) + TWO*POTCOEF(3,PPID)*R1 + &
          THREE*POTCOEF(4,PPID)*R1SQ + &
          FOUR*POTCOEF(5,PPID)*R1*R1SQ 
     !- POTCOEF(6,PPID)/R1SQ

     POTCOEF(12,PPID) = DPOLY*SCL_R1 + POTCOEF(7,PPID)*EXPTMP !+ &
     ! SIX*POTCOEF(8,PPID)/(R1*R6)

     DDPOLY = TWO*POTCOEF(3,PPID) + SIX*POTCOEF(4,PPID)*R1 + &
          TWELVE*POTCOEF(5,PPID)*R1SQ 
     !+ TWO*POTCOEF(6,PPID)/(R1*R1SQ)

     POTCOEF(13,PPID) = HALF*( (DPOLY*DPOLY + DDPOLY)*SCL_R1 + &
          POTCOEF(7,PPID)*POTCOEF(7,PPID)*EXPTMP) ! - &
     !SIX*SEVEN*POTCOEF(8,PPID)/(R1SQ*R6) )

     ! At the end of the join function:

     R1 = POTCOEF(9,PPID)
     DELTA = RCUT - R1
     DELTA2 = DELTA*DELTA
     DELTA3 = DELTA2*DELTA
     DELTA4 = DELTA3*DELTA

     POTCOEF(14,PPID) = (MINUSONE/DELTA3)*(THREE*POTCOEF(13,PPID)*DELTA2 + &
          SIX*POTCOEF(12,PPID)*DELTA + TEN*POTCOEF(11,PPID))

     POTCOEF(15,PPID) = (ONE/DELTA4)*(THREE*POTCOEF(13,PPID)*DELTA2 + &
          EIGHT*POTCOEF(12,PPID)*DELTA + FIFTEEN*POTCOEF(11,PPID))

     POTCOEF(16,PPID) = (MINUSONE/(TEN*DELTA3))*(SIX*POTCOEF(15,PPID)*DELTA2 + &
          THREE*POTCOEF(14,PPID)*DELTA + POTCOEF(13,PPID))

  ENDDO

  RETURN

END SUBROUTINE VDWTAILCOEF

