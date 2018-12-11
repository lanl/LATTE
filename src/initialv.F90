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

SUBROUTINE INITIALV

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K, N, CLOCK, COUNT
  REAL(LATTEPREC) :: RN(3)
  REAL(LATTEPREC) :: MYMASS, VELFACTOR
  REAL(LATTEPREC) :: CONV, TCURRENT, RESCALE
  REAL(LATTEPREC) :: MOMENTUM(3), STDDEV, BOLTZ
  REAL(LATTEPREC), EXTERNAL :: GAUSSRN
  IF (EXISTERROR) RETURN

  MOMENTUM = ZERO

  IF ( RNDIST .EQ. "UNIFORM" ) THEN

     DO I = 1, NATS

        MYMASS = MASS(ELEMPOINTER(I))

        VELFACTOR = SQRT(ONE/MYMASS)

        CALL RANDOM_NUMBER(RN)

        V(1,I) = (TWO*RN(1) - ONE) * VELFACTOR
        V(2,I) = (TWO*RN(2) - ONE) * VELFACTOR
        V(3,I) = (TWO*RN(3) - ONE) * VELFACTOR

        MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MYMASS
        MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MYMASS
        MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MYMASS

     ENDDO

  ELSEIF ( RNDIST .EQ. "GAUSSIAN" ) THEN

     SETTH = 0

     BOLTZ = ONE/KE2T

     DO I = 1, NATS

        MYMASS = MASS(ELEMPOINTER(I))
        STDDEV = SQRT(BOLTZ*TTARGET/(MYMASS*MVV2KE))

        V(1,I) = GAUSSRN(ZERO, STDDEV)
        V(2,I) = GAUSSRN(ZERO, STDDEV)
        V(3,I) = GAUSSRN(ZERO, STDDEV)

        MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MYMASS
        MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MYMASS
        MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MYMASS

     ENDDO

  ELSE

     CALL ERRORS("initialv","Choose either UNIFORM or GAUSSIAN for the &
          & random number distribution")

  ENDIF

  MOMENTUM = MOMENTUM/SUMMASS

  DO I = 1, NATS

     V(1,I) = V(1,I) - MOMENTUM(1)
     V(2,I) = V(2,I) - MOMENTUM(2)
     V(3,I) = V(3,I) - MOMENTUM(3)

  ENDDO

  CALL GETKE

  RESCALE = SQRT(TTARGET/TEMPERATURE)

  V = RESCALE * V

  RETURN

END SUBROUTINE INITIALV
