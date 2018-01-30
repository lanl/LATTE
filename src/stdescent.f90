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

SUBROUTINE STDESCENT(ITER, DELTAENERGY, DELTAF)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE RELAXCOMMON
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, ITER
  REAL(LATTEPREC) :: MAXF, DISP, MYSHIFT, DELTAENERGY, DELTAF
  REAL(LATTEPREC), PARAMETER :: MAXMCONST = 0.05D0, MINMCONST = 1.0D-3
  REAL(LATTEPREC), PARAMETER :: MAXSHIFT = 0.05D0
  IF (EXISTERROR) RETURN

  IF (ITER .EQ. 1) THEN

     RELCONST = TWO*MINMCONST

  ELSE

     ! If the energy and maximum force have gone down, increase the step size

     IF (DELTAENERGY .LE. ZERO .AND. DELTAF .LT. ZERO) THEN

        RELCONST = 1.01*RELCONST

     ELSEIF (DELTAENERGY .GT. ZERO .OR. DELTAF .GT. ZERO) THEN

        ! If the energy is increasing - decrease the step size

        RELCONST = 0.5*RELCONST

     ENDIF

  ENDIF

  RELCONST = MIN(RELCONST, MAXMCONST)
  RELCONST = MAX(RELCONST, MINMCONST)

  !  PREVF = MAXF

  DO I = 1, NATS

     DO J = 1, 3

        MYSHIFT = RELCONST*FTOT(J,I)

        IF (ABS(MYSHIFT) .GT. MAXSHIFT) MYSHIFT = SIGN(MAXSHIFT, MYSHIFT)

        CR(J,I) = CR(J,I) + MYSHIFT

     ENDDO

  ENDDO

  RETURN

END SUBROUTINE STDESCENT
