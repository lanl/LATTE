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

SUBROUTINE COULTAILCOEF

  USE CONSTANTS_MOD
  USE COULOMBARRAY
  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: DELTA, DELTA2, DELTA3, DELTA4
  REAL(LATTEPREC) :: R1, R12, R13
  IF (EXISTERROR) RETURN

  !
  ! The cut-off function looks like this:
  !
  ! t(R) = B1 + B2*(R - R1) + B3*(R - R1)^2 + B4*(R - R1)^3 + 
  !          B5*(R - R1)^4 + B6*(R - R1)^5
  !

  R1 = COULR1
  R12 = R1*COULR1
  R13 = R12*COULR1

  COULB(1) = ONE/R1

  COULB(2) = MINUSONE/R12

  COULB(3) = ONE/R13

  DELTA = COULCUT - COULR1
  DELTA2 = DELTA*DELTA
  DELTA3 = DELTA2*DELTA
  DELTA4 = DELTA3*DELTA

  COULB(4) = (MINUSONE/DELTA3)*(THREE*COULB(3)*DELTA2 + &
       SIX*COULB(2)*DELTA + TEN*COULB(1))

  COULB(5) = (ONE/DELTA4)*(THREE*COULB(3)*DELTA2 + &
       EIGHT*COULB(2)*DELTA + FIFTEEN*COULB(1))

  COULB(6) = (MINUSONE/(TEN*DELTA3))*(SIX*COULB(5)*DELTA2 + &
       THREE*COULB(4)*DELTA + COULB(3))

  RETURN

END SUBROUTINE COULTAILCOEF

