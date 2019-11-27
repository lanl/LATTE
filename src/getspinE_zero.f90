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

SUBROUTINE GETSPINEZERO

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY

  IMPLICIT NONE

  INTEGER :: I, J, K, MYINDEX
  REAL(LATTEPREC) :: DSPINS, DSPINP, DSPIND, DSPINF

  ESPIN_ZERO = ZERO
  MYINDEX = 0

  DO I = 1, NATS

     SELECT CASE(BASIS(ELEMPOINTER(I)))

     CASE("s")

        MYINDEX = MYINDEX + 1
        DSPINS = RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)
        ESPIN_ZERO = ESPIN_ZERO + WSS(ELEMPOINTER(I))*DSPINS*DSPINS

     CASE("p")

        MYINDEX = MYINDEX + 1
        DSPINS = RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)

        IF (BASIS(ELEMPOINTER(I)) .EQ. "s") THEN

           MYINDEX = MYINDEX + 1
           DSPINS = RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)

           ESPIN_ZERO = ESPIN_ZERO + WSS(ELEMPOINTER(I))*DSPINS*DSPINS

        ELSEIF (BASIS(ELEMPOINTER(I)) .EQ. "sp") THEN

           MYINDEX = MYINDEX + 1

           DSPINS = RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)

           DSPINP = ZERO

           MYINDEX = MYINDEX + 1
           DSPINP = DSPINP + RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)
           MYINDEX = MYINDEX + 1
           DSPINP = DSPINP + RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)
           MYINDEX = MYINDEX + 1
           DSPINP = DSPINP + RHOUPZERO(MYINDEX) - RHODOWNZERO(MYINDEX)

           ESPIN_ZERO = ESPIN_ZERO + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
                                !             TWO*WSP(ELEMPOINTER(I))*DSPINS*DSPINP + &
                WPP(ELEMPOINTER(I))*DSPINP*DSPINP

        ENDIF

     ENDDO

     ESPIN_ZERO = HALF*ESPIN_ZERO

     !  PRINT*, "ESPIN_ZERO= ", ESPIN_ZERO, DSPINP
     RETURN

   END SUBROUTINE GETSPINEZERO
