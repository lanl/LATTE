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

SUBROUTINE PARAFILEOPEN

#ifdef MPI_ON
  USE MPI
#endif

  USE CONSTANTS_MOD, ONLY: EXISTERROR

  IMPLICIT NONE

  INTEGER :: MYID, IERR, MYUNIT
  CHARACTER(LEN=50) :: FLNM

  IF (EXISTERROR) RETURN

#ifdef MPI_ON
  CALL MPI_COMM_RANK( MPI_COMM_WORLD, MYID, IERR )
#endif

  MYUNIT= 100 + MYID

  IF (MYID .LT. 10) THEN
     WRITE(FLNM,'(I1,"/pararep.dat")') MYID
  ELSEIF (MYID .GE. 10 .AND. MYID .LT. 100) THEN
     WRITE(FLNM,'(I2,"/pararep.dat")') MYID
  ELSEIF (MYID .GE. 100 .AND. MYID .LT. 1000) THEN
     WRITE(FLNM,'(I3,"/pararep.dat")') MYID
  ENDIF

  OPEN(UNIT=MYUNIT, STATUS="NEW", FILE=FLNM)

  RETURN

END SUBROUTINE PARAFILEOPEN

