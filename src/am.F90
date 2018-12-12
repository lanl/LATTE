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

FUNCTION AM(M, ALPHA)

  ! Build Am function defined in eqn. (5) of PRB 72 165107

  USE MYPRECISION

  IMPLICIT NONE

  REAL(LATTEPREC) :: AM, ALPHA
  !  REAL(LATTEPREC), EXTERNAL:: HEAVI
  INTEGER :: M

  ! Removed the call to the Heaviside step function to save on 
  ! cos or sin calculation in M .NE. 0

  ! MJC 1/10/14

  IF (M == 0) THEN
     AM = ONE / SQRT2
  ELSE IF ( M .GT. 0) THEN
     AM = REAL((-1)**M, LATTEPREC) * COS(ABS(M)*ALPHA)
  ELSE IF ( M .LT. 0) THEN
     AM = REAL((-1)**M, LATTEPREC) * SIN(ABS(M)*ALPHA)
  ENDIF

  RETURN

END FUNCTION AM
