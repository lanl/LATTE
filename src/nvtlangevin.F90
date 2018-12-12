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

!!$The algorithm developed by Gronbech-Jensen and Farago is implemented
!!$Molecular Physics 111 (2013) 983-991
!!$
!!$Contribution from Enrique Martinez

SUBROUTINE NVTLANGEVIN(ITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, N
  INTEGER :: ITER
  REAL(LATTEPREC) :: BETA, b, A
  REAL(LATTEPREC), EXTERNAL :: GAUSSRN
  REAL(LATTEPREC) :: MEAN, STDDEV
  REAL(LATTEPREC) :: BOLTZ, SQMVV2KE
  REAL(LATTEPREC) :: PREFRIC, PREB, V1TERM, V2TERM, V3TERM
  REAL(LATTEPREC) :: GAMMA1, GAMMA2, TSQRT, DTF, DTV, DTFM
  REAL(LATTEPREC) :: PREF1, PREF2, MELPREF, MULTFACT, MYMASS
  INTEGER :: OPTION
  IF (EXISTERROR) RETURN

  SETTH = 0

  BOLTZ = ONE/KE2T

  MEAN = ZERO

  DO I = 1, NATS

!!$     Update positions

     MYMASS = MASS(ELEMPOINTER(I))

     GAMMA1 = MYMASS/FRICTION

     STDDEV = SQRT(TWO*GAMMA1*BOLTZ*TTARGET*F2V*DT) 

     PREFRIC = (GAMMA1*DT)/(TWO*MYMASS)
     b = ONE / (ONE + PREFRIC)
     A = (ONE - PREFRIC)/(ONE + PREFRIC)
     PREB = (b*DT)/(TWO*MYMASS)

     CR(1,I) = CR(1,I) + b*DT*V(1,I) + PREB*DT*F2V*FTOT(1,I) + &
          PREB*GAUSSRN(MEAN,STDDEV)

     CR(2,I) = CR(2,I) + b*DT*V(2,I) + PREB*DT*F2V*FTOT(2,I) + &
          PREB*GAUSSRN(MEAN,STDDEV)

     CR(3,I) = CR(3,I) + b*DT*V(3,I) + PREB*DT*F2V*FTOT(3,I) + &
          PREB*GAUSSRN(MEAN,STDDEV)

     ! Update velocities

     MULTFACT = (DT*F2V)/(TWO*MYMASS)

     V(1,I) = A * (V(1,I) + MULTFACT * FTOT(1,I)) + &
          (b/MYMASS)*GAUSSRN(MEAN,STDDEV)

     V(2,I) = A * (V(2,I) + MULTFACT * FTOT(2,I)) + &
          (b/MYMASS)*GAUSSRN(MEAN,STDDEV)

     V(3,I) = A * (V(3,I) + MULTFACT * FTOT(3,I)) + &
          (b/MYMASS)*GAUSSRN(MEAN,STDDEV)

  ENDDO

  !
  ! Get new force to complete advance in V
  !

  CALL GETMDF(1, 100)

  DO I = 1, NATS

     MULTFACT = (DT*F2V)/(TWO*MASS(ELEMPOINTER(I)))

     V(1,I) = V(1,I) + MULTFACT*FTOT(1,I)
     V(2,I) = V(2,I) + MULTFACT*FTOT(2,I)
     V(3,I) = V(3,I) + MULTFACT*FTOT(3,I)

  ENDDO

  RETURN

END SUBROUTINE NVTLANGEVIN



