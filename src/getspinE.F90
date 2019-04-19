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

SUBROUTINE GETSPINE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE SPINARRAY

  IMPLICIT NONE

  INTEGER :: I, INDEX
  REAL(LATTEPREC) :: DSPINS, DSPINP, DSPIND, DSPINF
  IF (EXISTERROR) RETURN

  ESPIN = ZERO
  INDEX = 0

  DO I = 1, NATS

     SELECT CASE(BASIS(ELEMPOINTER(I)))

     CASE("s")

        DSPINS = DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS

     CASE("p")

        DSPINP = DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        ESPIN = ESPIN + WPP(ELEMPOINTER(I))*DSPINP*DSPINP

     CASE("d")

        DSPIND = DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        ESPIN = ESPIN + WDD(ELEMPOINTER(I))*DSPIND*DSPIND

     CASE("f")

        DSPINF = DELTASPIN(INDEX + 1)
        INDEX = INDEX + 1

        ESPIN = ESPIN + WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("sp")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPINP = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WPP(ELEMPOINTER(I))*DSPINP*DSPINP

     CASE("sd")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPIND = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND

     CASE("sf")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPINF = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF           

     CASE("pd")

        DSPINP = DELTASPIN(INDEX + 1)
        DSPIND = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND

     CASE("pf")

        DSPINP = DELTASPIN(INDEX + 1)
        DSPINF = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("df")

        DSPIND = DELTASPIN(INDEX + 1)
        DSPINF = DELTASPIN(INDEX + 2)
        INDEX = INDEX + 2

        ESPIN = ESPIN + WDD(ELEMPOINTER(I))*DSPIND*DSPIND + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("spd")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPINP = DELTASPIN(INDEX + 2)
        DSPIND = DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND

     CASE("spf")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPINP = DELTASPIN(INDEX + 2)
        DSPINF = DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("sdf")

        DSPINS = DELTASPIN(INDEX + 1)
        DSPIND = DELTASPIN(INDEX + 2)
        DSPINF = DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("pdf")

        DSPINP = DELTASPIN(INDEX + 1)
        DSPIND = DELTASPIN(INDEX + 2)
        DSPINF = DELTASPIN(INDEX + 3)
        INDEX = INDEX + 3

        ESPIN = ESPIN + WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     CASE("spdf") 

        DSPINS = DELTASPIN(INDEX + 1)
        DSPINP = DELTASPIN(INDEX + 2)
        DSPIND = DELTASPIN(INDEX + 3)
        DSPINF = DELTASPIN(INDEX + 4)
        INDEX = INDEX + 4

        ESPIN = ESPIN + WSS(ELEMPOINTER(I))*DSPINS*DSPINS + &
             WPP(ELEMPOINTER(I))*DSPINP*DSPINP + &
             WDD(ELEMPOINTER(I))*DSPIND*DSPIND + &
             WFF(ELEMPOINTER(I))*DSPINF*DSPINF

     END SELECT

  ENDDO

  ESPIN = HALF*ESPIN

  RETURN

END SUBROUTINE GETSPINE
