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

SUBROUTINE DUNIVSCALE_SUB(R, A, DC, X, DXDR)

  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: A(14)
  REAL(LATTEPREC) :: DC(3), DXDR(3), TMPDXDR, RMOD
  REAL(LATTEPREC) :: X, R, RMINUSR1, DPOLYNOM

  CALL UNIVSCALE_SUB(R, A, X)

  IF (R .LE. A(7)) THEN

     RMOD = R - A(6)
     DPOLYNOM = A(2) + RMOD*(TWO*A(3) + RMOD*(THREE*A(4) + FOUR*A(5)*RMOD)) 

     TMPDXDR = DPOLYNOM*X

  ELSEIF (R .GT. A(7) .AND. R .LT. A(8)) THEN

     RMINUSR1 = R - A(7)

     TMPDXDR = A(10) + RMINUSR1*(TWO*A(11) + &
          RMINUSR1*(THREE*A(12) + RMINUSR1*(FOUR*A(13) + &
          RMINUSR1*FIVE*A(14))))

     TMPDXDR = A(1)*TMPDXDR

  ELSE

     TMPDXDR = ZERO

  END IF

  DXDR = -TMPDXDR*DC

  RETURN

END SUBROUTINE DUNIVSCALE_SUB
