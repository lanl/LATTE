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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\

!!$ Nose thermostat from Anders Niklasson derivation
!!$ Contributing author: Enrique Martinez

SUBROUTINE NVTNH

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, N
  INTEGER :: ITER
  REAL(LATTEPREC) :: KEV, KO
  REAL(LATTEPREC) :: PREF1, PREF2, NOSEPREF
  IF (EXISTERROR) RETURN

  KEV = ZERO
  DO I = 1, NATS

     KEV = KEV + MASS(ELEMPOINTER(I)) * & 
         & (V(1,I)*V(1,I) + V(2,I)*V(2,I) + V(3,I)*V(3,I))

  ENDDO

  KEV = MVV2KE*KEV/TWO

  KO = THREE*REAL(NATS)*TTARGET/(TWO*KE2T)

  DGAMMA = DGAMMA + (DT/FRICTION)*(KEV - KO)/TWO
  GAMMA = GAMMA + DT*DGAMMA/TWO
  DGAMMA = DGAMMA + (DT/FRICTION)*(KEV - KO)/TWO

  PREF1 = F2V*DT/TWO
  NOSEPREF = ONE / (ONE + (DT*DGAMMA/TWO))

  DO I = 1, NATS

     PREF2 = PREF1/MASS(ELEMPOINTER(I))

     V(1,I) = NOSEPREF*( V(1,I) + PREF2*FTOT(1,I) )
     V(2,I) = NOSEPREF*( V(2,I) + PREF2*FTOT(2,I) )
     V(3,I) = NOSEPREF*( V(3,I) + PREF2*FTOT(3,I) )

     CR(1,I) = CR(1,I) + DT*V(1,I)
     CR(2,I) = CR(2,I) + DT*V(2,I)
     CR(3,I) = CR(3,I) + DT*V(3,I)

  ENDDO

  CALL GETMDF(1, 100)


  NOSEPREF = ONE - DT*DGAMMA/TWO

  DO I = 1, NATS

     PREF2 = PREF1/MASS(ELEMPOINTER(I))

     V(1,I) = NOSEPREF*V(1,I) + PREF2*FTOT(1,I)
     V(2,I) = NOSEPREF*V(2,I) + PREF2*FTOT(2,I)
     V(3,I) = NOSEPREF*V(3,I) + PREF2*FTOT(3,I)

  ENDDO

  KEV = ZERO
  DO I = 1, NATS
     KEV = KEV + MASS(ELEMPOINTER(I)) * & 
        &  (V(1,I)*V(1,I) + V(2,I)*V(2,I) + V(3,I)*V(3,I))      
  ENDDO

  KEV = MVV2KE*KEV/TWO

  DGAMMA = DGAMMA + (DT/FRICTION)*(KEV - KO)/TWO
  GAMMA = GAMMA + DT*DGAMMA/TWO
  DGAMMA = DGAMMA + (DT/FRICTION)*(KEV - KO)/TWO

  CONSMOT = FRICTION*DGAMMA*DGAMMA/TWO + TWO*KO*GAMMA

  RETURN

END SUBROUTINE NVTNH
