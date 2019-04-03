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

SUBROUTINE UNIVSCALE_SUB(R, A, X)

  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: A(14), X, R, RMINUSR1, POLYNOM, RMOD

  IF (R .LE. A(7)) THEN

     RMOD = R - A(6)

     POLYNOM = RMOD*(A(2) + RMOD*(A(3) + RMOD*(A(4) + A(5)*RMOD)))

     X = EXP(POLYNOM)

  ELSEIF (R .GT. A(7) .AND. R .LT. A(8)) THEN

     RMINUSR1 = R - A(7)

     X = A(9) + RMINUSR1*(A(10) + &
          RMINUSR1*(A(11) + RMINUSR1*(A(12) + &
          RMINUSR1*(A(13) + RMINUSR1*A(14)))))

  ELSE

     X = ZERO

  END IF

  X = A(1)*X

  RETURN

END SUBROUTINE UNIVSCALE_SUB
