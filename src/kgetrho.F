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

SUBROUTINE KGETRHO

  USE CONSTANTS_MOD
  USE MYPRECISION
  USE KSPACEARRAY
  USE DIAGARRAY
  !  USE MPI

  IMPLICIT NONE
  IF (EXISTERROR) RETURN

  !  INTEGER :: START_CLOCK, STOP_CLOCK, CLOCK_RATE, CLOCK_MAX

  !  INTEGER :: I, J, MYID, IERR
  ! Lets first build our NKTOT Hamiltonians (we'll need all of them
  ! in order to get the density matrix via SP2....

  !  CALL MPI_COMM_RANK(MPI_COMM_WORLD, MYID, IERR)

  !  CALL SYSTEM_CLOCK(START_CLOCK, CLOCK_RATE, CLOCK_MAX)

  CALL KDIAGMYH

  !  CALL SYSTEM_CLOCK(STOP_CLOCK, CLOCK_RATE, CLOCK_MAX)

  !  PRINT*, "# Diag time (s) = ", &
  !      REAL(STOP_CLOCK - START_CLOCK)/REAL(CLOCK_RATE)

  !  DO I = 1, NKTOT
  !     PRINT*, I, MYID, KEVALS(1,I)
  !  ENDDO

  CALL KBOEVECS



  !  DO I = 1, NKTOT
  !        PRINT*, I, MYID, KBO(1,1,I), KBO(HDIM, HDIM, I)
  !  ENDDO

  RETURN

END SUBROUTINE KGETRHO
