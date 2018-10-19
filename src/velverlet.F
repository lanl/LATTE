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

SUBROUTINE VELVERLET(CURRITER)

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MDARRAY
  USE MYPRECISION
  USE CONSTRAINTS_MOD

  IMPLICIT NONE

  INTEGER :: I, J, K, CURRITER
  REAL(LATTEPREC) :: PREFACTOR, PREPREFACT 
  IF (EXISTERROR) RETURN

  PREPREFACT = HALF*F2V*DT

  !
  ! Half timestep advance in V
  !

  !  IF (FREEZE .EQ. 1) CALL FREEZE_ATOMS(V,FTOT)

  DO I = 1, NATS

     PREFACTOR = PREPREFACT/MASS(ELEMPOINTER(I))

     V(1,I) = V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO

  CR = CR + DT*V

  !
  ! Get new force to complete advance in V
  !

  CALL GETMDF(1, CURRITER)

  !
  ! Now finish advancing V with F(t + dt)
  !

  DO I = 1, NATS

     PREFACTOR = PREPREFACT/MASS(ELEMPOINTER(I))

     V(1,I) = V(1,I) + PREFACTOR*FTOT(1,I)
     V(2,I) = V(2,I) + PREFACTOR*FTOT(2,I)
     V(3,I) = V(3,I) + PREFACTOR*FTOT(3,I)

  ENDDO

  RETURN

END SUBROUTINE VELVERLET

