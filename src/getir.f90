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

SUBROUTINE GETIR

  USE CONSTANTS_MOD
  USE SETUPARRAY
  USE MYPRECISION
  USE MDARRAY

  IMPLICIT NONE
  IF (EXISTERROR) RETURN


  FTOT = ZERO

  !Change in coordinate x of atom 1
  CR(1,1) = CR(1,1) + 0.001d0

  !Change in coordinate y of atom 1
  CR(2,1) = CR(2,1) + 0.001d0
  
  !Computing forces
  CALL GETFORCE

  WRITE(*,*)"Force of atom 1 in x", FTOT(1,1)
  WRITE(*,*)"Force of atom 1 in y", FTOT(2,1)


  WRITE(*,*)"Mass of atom 1",  MASS(ELEMPOINTER(1))

  WRITE(*,*)"==========================="
  WRITE(*,*)"Im stoping in subroutine getir"
  WRITE(*,*)"Grep for GETIR"
  WRITE(*,*)"==========================="
  
  STOP

  RETURN

END SUBROUTINE GETIR
