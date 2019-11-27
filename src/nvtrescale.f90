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

SUBROUTINE NVTRESCALE

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION

  IMPLICIT NONE

  INTEGER :: I, J, K
  REAL(LATTEPREC) :: PREFACTOR, PREPREFACT
  REAL(LATTEPREC) :: CHI, TMPTEMP, MYTEMP, MOMENTUM(3)
  REAL(LATTEPREC) :: MAXCHI = 1.05D0
  IF (EXISTERROR) RETURN

  PREPREFACT = HALF*F2V*DT

  TMPTEMP = ZERO

  DO I = 1, AVEPER
     TMPTEMP = TMPTEMP + THIST(I)
  ENDDO

  MYTEMP = TMPTEMP/FLOAT(AVEPER)

  !  CHI = SQRT(TTARGET/MYTEMP)

  ! CHI = MIN(CHI, MAXCHI)
  CHI = MIN(SQRT(TTARGET/MYTEMP), MAXCHI)

  !  PRINT*, "MYTEMP = ", MYTEMP, "CHI = ", CHI

  DO I = 1, NATS

     PREFACTOR = PREPREFACT/MASS(ELEMPOINTER(I))

     !
     ! Half timestep advance in V with velocity rescale
     !

     V(1,I) = CHI*V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = CHI*V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = CHI*V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO

  !
  ! Whole timestep advance in positions
  !

  CR = CR + DT*V

  !
  ! Get new force to complete advance in V
  !

  CALL GETMDF(1, 100)

  !
  ! Now finish advancing V with F(t + dt)
  !

  DO I = 1, NATS

     PREFACTOR = PREPREFACT/MASS(ELEMPOINTER(I))

     V(1,I) = V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO

  ! Remove drift

  MOMENTUM = ZERO

  DO I = 1, NATS
     MOMENTUM(1) = MOMENTUM(1) + V(1,I)*MASS(ELEMPOINTER(I))
     MOMENTUM(2) = MOMENTUM(2) + V(2,I)*MASS(ELEMPOINTER(I))
     MOMENTUM(3) = MOMENTUM(3) + V(3,I)*MASS(ELEMPOINTER(I))
  ENDDO

  MOMENTUM = MOMENTUM/SUMMASS

  DO I = 1, NATS
     V(1,I) = V(1,I) - MOMENTUM(1)
     V(2,I) = V(2,I) - MOMENTUM(2)
     V(3,I) = V(3,I) - MOMENTUM(3)
  ENDDO

  RETURN

END SUBROUTINE NVTRESCALE

